#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "HeterogeneousCore/SonicTriton/interface/TritonClient.h"
#include "HeterogeneousCore/SonicTriton/interface/triton_utils.h"

#include "grpc_client.h"
#include "grpc_service.pb.h"

#include <string>
#include <cmath>
#include <chrono>
#include <exception>
#include <sstream>
#include <utility>
#include <tuple>

namespace ni = nvidia::inferenceserver;
namespace nic = ni::client;

//based on https://github.com/triton-inference-server/server/blob/v2.3.0/src/clients/c++/examples/simple_grpc_async_infer_client.cc
//and https://github.com/triton-inference-server/server/blob/v2.3.0/src/clients/c++/perf_client/perf_client.cc

TritonClient::TritonClient(const edm::ParameterSet& params)
    : SonicClient(params),
      verbose_(params.getUntrackedParameter<bool>("verbose")),
      options_(params.getParameter<std::string>("modelName")) {
  clientName_ = "TritonClient";
  //will get overwritten later, just used in constructor
  fullDebugName_ = clientName_;

  //connect to the server
  //TODO: add SSL options
  std::string url(params.getUntrackedParameter<std::string>("address") + ":" +
                  std::to_string(params.getUntrackedParameter<unsigned>("port")));
  triton_utils::throwIfError(nic::InferenceServerGrpcClient::Create(&client_, url, false),
                             "TritonClient(): unable to create inference context");

  //set options
  options_.model_version_ = params.getParameter<std::string>("modelVersion");
  //convert seconds to microseconds
  options_.client_timeout_ = params.getUntrackedParameter<unsigned>("timeout") * 1e6;

  //config needed for batch size
  inference::ModelConfigResponse modelConfigResponse;
  triton_utils::throwIfError(client_->ModelConfig(&modelConfigResponse, options_.model_name_, options_.model_version_),
                             "TritonClient(): unable to get model config");
  inference::ModelConfig modelConfig(modelConfigResponse.config());

  //check batch size limitations (after i/o setup)
  //triton uses max batch size = 0 to denote a model that does not support batching
  //but for models that do support batching, a given event may set batch size 0 to indicate no valid input is present
  //so set the local max to 1 and keep track of "no batch" case
  maxBatchSize_ = modelConfig.max_batch_size();
  noBatch_ = maxBatchSize_ == 0;
  maxBatchSize_ = std::max(1u, maxBatchSize_);

  //get model info
  inference::ModelMetadataResponse modelMetadata;
  triton_utils::throwIfError(client_->ModelMetadata(&modelMetadata, options_.model_name_, options_.model_version_),
                             "TritonClient(): unable to get model metadata");

  //get input and output (which know their sizes)
  const auto& nicInputs = modelMetadata.inputs();
  const auto& nicOutputs = modelMetadata.outputs();

  //report all model errors at once
  std::stringstream msg;
  std::string msg_str;

  //currently no use case is foreseen for a model with zero inputs or outputs
  if (nicInputs.empty())
    msg << "Model on server appears malformed (zero inputs)\n";

  if (nicOutputs.empty())
    msg << "Model on server appears malformed (zero outputs)\n";

  //stop if errors
  msg_str = msg.str();
  if (!msg_str.empty())
    throw cms::Exception("ModelErrors") << msg_str;

  //setup input map
  std::stringstream io_msg;
  if (verbose_)
    io_msg << "Model inputs: "
           << "\n";
  inputsTriton_.reserve(nicInputs.size());
  for (const auto& nicInput : nicInputs) {
    const auto& iname = nicInput.name();
    auto [curr_itr, success] = input_.emplace(
        std::piecewise_construct, std::forward_as_tuple(iname), std::forward_as_tuple(iname, nicInput, noBatch_));
    auto& curr_input = curr_itr->second;
    inputsTriton_.push_back(curr_input.data());
    if (verbose_) {
      io_msg << "  " << iname << " (" << curr_input.dname() << ", " << curr_input.byteSize()
             << " b) : " << triton_utils::printColl(curr_input.shape()) << "\n";
    }
  }

  //allow selecting only some outputs from server
  const auto& v_outputs = params.getUntrackedParameter<std::vector<std::string>>("outputs");
  std::unordered_set s_outputs(v_outputs.begin(), v_outputs.end());

  //setup output map
  if (verbose_)
    io_msg << "Model outputs: "
           << "\n";
  outputsTriton_.reserve(nicOutputs.size());
  for (const auto& nicOutput : nicOutputs) {
    const auto& oname = nicOutput.name();
    if (!s_outputs.empty() and s_outputs.find(oname) == s_outputs.end())
      continue;
    auto [curr_itr, success] = output_.emplace(
        std::piecewise_construct, std::forward_as_tuple(oname), std::forward_as_tuple(oname, nicOutput, noBatch_));
    auto& curr_output = curr_itr->second;
    outputsTriton_.push_back(curr_output.data());
    if (verbose_) {
      io_msg << "  " << oname << " (" << curr_output.dname() << ", " << curr_output.byteSize()
             << " b) : " << triton_utils::printColl(curr_output.shape()) << "\n";
    }
    if (!s_outputs.empty())
      s_outputs.erase(oname);
  }

  //check if any requested outputs were not available
  if (!s_outputs.empty())
    throw cms::Exception("MissingOutput")
        << "Some requested outputs were not available on the server: " << triton_utils::printColl(s_outputs);

  //check requested batch size and propagate to inputs and outputs
  setBatchSize(params.getUntrackedParameter<unsigned>("batchSize"));

  //print model info
  std::stringstream model_msg;
  if (verbose_) {
    model_msg << "Model name: " << options_.model_name_ << "\n"
              << "Model version: " << options_.model_version_ << "\n"
              << "Model max batch size: " << (noBatch_ ? 0 : maxBatchSize_) << "\n";
    edm::LogInfo(fullDebugName_) << model_msg.str() << io_msg.str();
  }
}

bool TritonClient::setBatchSize(unsigned bsize) {
  if (bsize > maxBatchSize_) {
    edm::LogWarning(fullDebugName_) << "Requested batch size " << bsize << " exceeds server-specified max batch size "
                                    << maxBatchSize_ << ". Batch size will remain as" << batchSize_;
    return false;
  } else {
    batchSize_ = bsize;
    //set for input and output
    for (auto& element : input_) {
      element.second.setBatchSize(bsize);
    }
    for (auto& element : output_) {
      element.second.setBatchSize(bsize);
    }
    return true;
  }
}

void TritonClient::reset() {
  for (auto& element : input_) {
    element.second.reset();
  }
  for (auto& element : output_) {
    element.second.reset();
  }
}

bool TritonClient::getResults(std::shared_ptr<nic::InferResult> results) {
  for (auto& [oname, output] : output_) {
    //set shape here before output becomes const
    if (output.variableDims()) {
      std::vector<int64_t> tmp_shape;
      bool status = triton_utils::warnIfError(results->Shape(oname, &tmp_shape),
                                              "getResults(): unable to get output shape for " + oname);
      if (!status)
        return status;
      output.setShape(tmp_shape, false);
    }
    //extend lifetime
    output.setResult(results);
  }

  return true;
}

//default case for sync and pseudo async
void TritonClient::evaluate() {
  //in case there is nothing to process
  if (batchSize_ == 0) {
    finish(true);
    return;
  }

  // Get the status of the server prior to the request being made.
  const auto& start_status = getServerSideStatus();

  if (mode_ == SonicMode::Async) {
    //non-blocking call
    auto t1 = std::chrono::high_resolution_clock::now();
    bool status = triton_utils::warnIfError(
        client_->AsyncInfer(
            [t1, start_status, this](nic::InferResult* results) {
              //get results
              std::shared_ptr<nic::InferResult> results_ptr(results);
              bool status = triton_utils::warnIfError(results_ptr->RequestStatus(), "evaluate(): unable to get result");
              if (!status) {
                finish(false);
                return;
              }
              auto t2 = std::chrono::high_resolution_clock::now();

              if (!debugName_.empty())
                edm::LogInfo(fullDebugName_)
                    << "Remote time: " << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

              const auto& end_status = getServerSideStatus();

              if (verbose()) {
                const auto& stats = summarizeServerStats(start_status, end_status);
                reportServerSideStats(stats);
              }

              //check result
              status = getResults(results_ptr);

              //finish
              finish(status);
            },
            options_,
            inputsTriton_,
            outputsTriton_),
        "evaluate(): unable to launch async run");

    //if AsyncRun failed, finish() wasn't called
    if (!status)
      finish(false);
  } else {
    //blocking call
    auto t1 = std::chrono::high_resolution_clock::now();
    nic::InferResult* results;
    bool status = triton_utils::warnIfError(client_->Infer(&results, options_, inputsTriton_, outputsTriton_),
                                            "evaluate(): unable to run and/or get result");
    if (!status) {
      finish(false);
      return;
    }

    auto t2 = std::chrono::high_resolution_clock::now();
    if (!debugName_.empty())
      edm::LogInfo(fullDebugName_) << "Remote time: "
                                   << std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

    const auto& end_status = getServerSideStatus();

    if (verbose()) {
      const auto& stats = summarizeServerStats(start_status, end_status);
      reportServerSideStats(stats);
    }

    std::shared_ptr<nic::InferResult> results_ptr(results);
    status = getResults(results_ptr);

    finish(status);
  }
}

void TritonClient::reportServerSideStats(const TritonClient::ServerSideStats& stats) const {
  std::stringstream msg;

  // https://github.com/triton-inference-server/server/blob/v2.3.0/src/clients/c++/perf_client/inference_profiler.cc
  const uint64_t count = stats.success_count_;
  msg << "  Inference count: " << stats.inference_count_ << "\n";
  msg << "  Execution count: " << stats.execution_count_ << "\n";
  msg << "  Successful request count: " << count << "\n";

  if (count > 0) {
    auto get_avg_us = [count](uint64_t tval) {
      constexpr uint64_t us_to_ns = 1000;
      return tval / us_to_ns / count;
    };

    const uint64_t cumm_avg_us = get_avg_us(stats.cumm_time_ns_);
    const uint64_t queue_avg_us = get_avg_us(stats.queue_time_ns_);
    const uint64_t compute_input_avg_us = get_avg_us(stats.compute_input_time_ns_);
    const uint64_t compute_infer_avg_us = get_avg_us(stats.compute_infer_time_ns_);
    const uint64_t compute_output_avg_us = get_avg_us(stats.compute_output_time_ns_);
    const uint64_t compute_avg_us = compute_input_avg_us + compute_infer_avg_us + compute_output_avg_us;
    const uint64_t overhead =
        (cumm_avg_us > queue_avg_us + compute_avg_us) ? (cumm_avg_us - queue_avg_us - compute_avg_us) : 0;

    msg << "  Avg request latency: " << cumm_avg_us << " usec"
        << "\n"
        << "  (overhead " << overhead << " usec + "
        << "queue " << queue_avg_us << " usec + "
        << "compute input " << compute_input_avg_us << " usec + "
        << "compute infer " << compute_infer_avg_us << " usec + "
        << "compute output " << compute_output_avg_us << " usec)" << std::endl;
  }

  if (!debugName_.empty())
    edm::LogInfo(fullDebugName_) << msg.str();
}

TritonClient::ServerSideStats TritonClient::summarizeServerStats(const inference::ModelStatistics& start_status,
                                                                 const inference::ModelStatistics& end_status) const {
  TritonClient::ServerSideStats server_stats;

  server_stats.inference_count_ = end_status.inference_count() - start_status.inference_count();
  server_stats.execution_count_ = end_status.execution_count() - start_status.execution_count();
  server_stats.success_count_ =
      end_status.inference_stats().success().count() - start_status.inference_stats().success().count();
  server_stats.cumm_time_ns_ =
      end_status.inference_stats().success().ns() - start_status.inference_stats().success().ns();
  server_stats.queue_time_ns_ = end_status.inference_stats().queue().ns() - start_status.inference_stats().queue().ns();
  server_stats.compute_input_time_ns_ =
      end_status.inference_stats().compute_input().ns() - start_status.inference_stats().compute_input().ns();
  server_stats.compute_infer_time_ns_ =
      end_status.inference_stats().compute_infer().ns() - start_status.inference_stats().compute_infer().ns();
  server_stats.compute_output_time_ns_ =
      end_status.inference_stats().compute_output().ns() - start_status.inference_stats().compute_output().ns();

  return server_stats;
}

inference::ModelStatistics TritonClient::getServerSideStatus() const {
  if (verbose_) {
    inference::ModelStatisticsResponse resp;
    triton_utils::warnIfError(client_->ModelInferenceStatistics(&resp, options_.model_name_, options_.model_version_),
                              "getServerSideStatus(): unable to get model statistics");
    return *(resp.model_stats().begin());
  }
  return inference::ModelStatistics{};
}

//for fillDescriptions
void TritonClient::fillPSetDescription(edm::ParameterSetDescription& iDesc) {
  edm::ParameterSetDescription descClient;
  fillBasePSetDescription(descClient);
  descClient.add<std::string>("modelName");
  descClient.add<std::string>("modelVersion", "");
  //server parameters should not affect the physics results
  descClient.addUntracked<unsigned>("batchSize");
  descClient.addUntracked<std::string>("address");
  descClient.addUntracked<unsigned>("port");
  descClient.addUntracked<unsigned>("timeout");
  descClient.addUntracked<bool>("verbose", false);
  descClient.addUntracked<std::vector<std::string>>("outputs", {});
  iDesc.add<edm::ParameterSetDescription>("Client", descClient);
}
