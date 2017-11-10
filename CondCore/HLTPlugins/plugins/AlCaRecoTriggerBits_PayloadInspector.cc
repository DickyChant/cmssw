#include "CondCore/Utilities/interface/PayloadInspectorModule.h"
#include "CondCore/Utilities/interface/PayloadInspector.h"
#include "CondCore/CondDB/interface/Time.h"
#include "CondFormats/HLTObjects/interface/AlCaRecoTriggerBits.h"

#include <memory>
#include <sstream>
#include <iostream>
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"

namespace {
  
  /************************************************
    Display AlCaRecoTriggerBits mapping
  *************************************************/
  class AlCaRecoTriggerBits_Display: public cond::payloadInspector::PlotImage<AlCaRecoTriggerBits> {
  public:
    AlCaRecoTriggerBits_Display() : cond::payloadInspector::PlotImage<AlCaRecoTriggerBits>( "Table of AlCaRecoTriggerBits" ){
      setSingleIov( true );
    }

    bool fill( const std::vector<std::tuple<cond::Time_t,cond::Hash> >& iovs ) override{
      auto iov = iovs.front();
      std::shared_ptr<AlCaRecoTriggerBits> payload = fetchPayload( std::get<1>(iov) );

      std::string IOVsince  = std::to_string(std::get<0>(iov));

      // Get map of strings to concatenated list of names of HLT paths:
      typedef std::map<std::string, std::string> TriggerMap;
      const TriggerMap &triggerMap = payload->m_alcarecoToTrig;

      unsigned int mapsize = triggerMap.size();
      float pitch = 1./(mapsize*1.1);

      float y, x1, x2;
      std::vector<float> y_x1,y_x2,y_line;
      std::vector<std::string>   s_x1,s_x2,s_x3;

      y = 1.0; x1 = 0.02; x2 = x1+0.30;
	
      y -= pitch; 
      y_x1.push_back(y);
      s_x1.push_back("#scale[1.2]{Key}");
      y_x2.push_back(y);
      s_x2.push_back("#scale[1.2]{in IOV: "+IOVsince+"}");

      y -= pitch/2.;
      y_line.push_back(y);

      for(const auto &element : triggerMap){
	//std::cout<< element.first << " : " ;

	y -= pitch; 
	y_x1.push_back(y);
	s_x1.push_back(element.first.c_str());

	std::map<int,std::string> output;
	int count=0;
	const std::vector<std::string> paths = payload->decompose(element.second);
	for (unsigned int iPath = 0; iPath < paths.size(); ++iPath) {
	  if((output[count]+paths[iPath]).length()<60){
	    output[count]+=paths[iPath];
	    output[count]+=";";
	  } else {
	    count++;
	    output[count]+=paths[iPath];
	  }
	}
      	
	for (unsigned int br=0; br<output.size();br++){
	  y_x2.push_back(y);
	  s_x2.push_back("#color[2]{"+output[br]+"}");
	  if(br!=output.size()-1) y-=pitch;
	}

	y_line.push_back(y-(pitch/2.));
       
	//std::cout << std::endl;
      }

      TCanvas canvas("AlCaRecoTriggerBits","AlCaRecoTriggerBits",2000,std::max(y_x1.size(),y_x2.size())*40); 
      TLatex l;
      // Draw the columns titles
      l.SetTextAlign(12);

      float newpitch = 1/(std::max(y_x1.size(),y_x2.size())*1.1);
      float factor  = newpitch/pitch;
      l.SetTextSize(newpitch-0.002);
      canvas.cd();
      for(unsigned int i=0;i<y_x1.size();i++){
	l.DrawLatexNDC(x1,1-(1-y_x1[i])*factor,s_x1[i].c_str());
      }

      for(unsigned int i=0;i<y_x2.size();i++){
	l.DrawLatexNDC(x2,1-(1-y_x2[i])*factor,s_x2[i].c_str());
      }  

      canvas.cd();
      canvas.Update();

      TLine lines[y_line.size()];
      unsigned int iL=0;
      for (const auto & line : y_line){
	//std::cout<<1-(1-line)*factor<<std::endl;
	lines[iL] = TLine(gPad->GetUxmin(),1-(1-line)*factor,gPad->GetUxmax(),1-(1-line)*factor);
	lines[iL].SetLineWidth(1);
	lines[iL].SetLineStyle(9);
	lines[iL].SetLineColor(2);
	lines[iL].Draw("same");
	iL++;
      }

      std::string fileName(m_imageFileName);
      canvas.SaveAs(fileName.c_str());
      return true;
    }
  };

  /************************************************
    Compare AlCaRecoTriggerBits mapping
  *************************************************/
  class AlCaRecoTriggerBits_Compare: public cond::payloadInspector::PlotImage<AlCaRecoTriggerBits> {
  public:
    AlCaRecoTriggerBits_Compare() : cond::payloadInspector::PlotImage<AlCaRecoTriggerBits>( "Table of AlCaRecoTriggerBits comparison" ){
      setSingleIov( false );
    }

    bool fill( const std::vector<std::tuple<cond::Time_t,cond::Hash> >& iovs ) override{

      std::vector<std::tuple<cond::Time_t,cond::Hash> > sorted_iovs = iovs;
      
      // make absolute sure the IOVs are sortd by since
      std::sort(begin(sorted_iovs), end(sorted_iovs), [](auto const &t1, auto const &t2) {
	  return std::get<0>(t1) < std::get<0>(t2);
	});
      
      auto firstiov  = sorted_iovs.front();
      auto lastiov   = sorted_iovs.back();
      
      std::shared_ptr<AlCaRecoTriggerBits> last_payload  = fetchPayload( std::get<1>(lastiov) );
      std::shared_ptr<AlCaRecoTriggerBits> first_payload = fetchPayload( std::get<1>(firstiov) );
      
      std::string lastIOVsince  = std::to_string(std::get<0>(lastiov));
      std::string firstIOVsince = std::to_string(std::get<0>(firstiov));
      
      // Get map of strings to concatenated list of names of HLT paths:
      typedef std::map<std::string, std::string> TriggerMap;
      const TriggerMap &first_triggerMap = first_payload->m_alcarecoToTrig;
      const TriggerMap &last_triggerMap  = last_payload->m_alcarecoToTrig;

      std::vector<std::string> first_keys, not_in_first_keys;
      std::vector<std::string> last_keys, not_in_last_keys;

      // fill the vector of first keys
      for (const auto& element : first_triggerMap){
	first_keys.push_back(element.first);
      }

      // fill the vector of last keys
      for (const auto& element : last_triggerMap){
	last_keys.push_back(element.first);
      }

      // find the elements not in common
      std::set_difference(first_keys.begin(),first_keys.end(),last_keys.begin(),last_keys.end(), 
			  std::inserter(not_in_last_keys, not_in_last_keys.begin()));
      
      std::set_difference(last_keys.begin(),last_keys.end(),first_keys.begin(),first_keys.end(), 
			  std::inserter(not_in_first_keys, not_in_first_keys.begin()));

      float pitch = 0.013;
      float y, x1, x2, x3;

      std::vector<float>         y_x1,y_x2,y_x3,y_line;
      std::vector<std::string>   s_x1,s_x2,s_x3;

      y  = 1.0; 
      x1 = 0.02; 
      x2 = x1+0.20; 
      x3 = x2+0.30;

      y -= pitch; 
      y_x1.push_back(y);
      s_x1.push_back("#scale[1.2]{Key}");
      y_x2.push_back(y);
      s_x2.push_back("#scale[1.2]{in IOV: "+firstIOVsince+"}");
      y_x3.push_back(y);
      s_x3.push_back("#scale[1.2]{in IOV: "+lastIOVsince+"}");
      y -= pitch/3;

      // print the ones missing in the last key
      for(const auto& key : not_in_last_keys ) {
	//std::cout<< key ;
	y -= pitch;  
	y_x1.push_back(y);
	s_x1.push_back(key);

	const std::vector<std::string> missing_in_last_paths = first_payload->decompose(first_triggerMap.at(key));
	
	std::map<int,std::string> output;
	int count=0;
	for (unsigned int iPath = 0; iPath < missing_in_last_paths.size(); ++iPath) {
	  //std::cout << missing_in_last_paths[iPath] << " ; " ;

	  if((output[count]+missing_in_last_paths[iPath]).length()<60){
	    output[count]+=missing_in_last_paths[iPath];
	    output[count]+=";";
	  } else {
	    count++;
	    output[count]+=missing_in_last_paths[iPath];
	  }
	}
	
	for (unsigned int br=0; br<output.size();br++){
	  y_x2.push_back(y);
	  s_x2.push_back("#color[2]{"+output[br]+"}");
	  if(br!=output.size()-1) y -=pitch;
	}
	//std::cout << " |||||| not in last";
	//std::cout << std::endl;

	y_line.push_back(y-0.008);

      }
      

      // print the ones missing in the first key
      for(const auto& key : not_in_first_keys ) {
	//std::cout<< key ;
	y -= pitch;  
	y_x1.push_back(y);
	s_x1.push_back(key);
	const std::vector<std::string> missing_in_first_paths = last_payload->decompose(last_triggerMap.at(key));

	//std::cout << " not in first ||||||";
	std::map<int,std::string> output;
	int count=0;
	for (unsigned int iPath = 0; iPath < missing_in_first_paths.size(); ++iPath) {
	  //std::cout << missing_in_first_paths[iPath] << " ; " ;

	  if((output[count]+missing_in_first_paths[iPath]).length()<60){
	    output[count]+=missing_in_first_paths[iPath];
	    output[count]+=";";
	  } else {
	    count++;
	    output[count]+=missing_in_first_paths[iPath];
	  }
	}

	for (unsigned int br=0; br<output.size();br++){
	  y_x3.push_back(y);
	  s_x3.push_back("#color[4]{"+output[br]+"}");
	  if(br!=output.size()-1) y -= pitch;
	}
	//std::cout << std::endl;
	
	y_line.push_back(y-0.008);
	
      }

      for(const auto &element : first_triggerMap){

      	if(last_triggerMap.find(element.first)!=last_triggerMap.end()){

      	  auto lastElement = last_triggerMap.find(element.first);
	
      	  std::string output;
      	  const std::vector<std::string> first_paths = first_payload->decompose(element.second);
      	  const std::vector<std::string> last_paths  = last_payload->decompose(lastElement->second);

      	  std::vector<std::string> not_in_first;
      	  std::vector<std::string> not_in_last;

      	  std::set_difference(first_paths.begin(),first_paths.end(),last_paths.begin(),last_paths.end(), 
      			      std::inserter(not_in_last, not_in_last.begin()));

      	  std::set_difference(last_paths.begin(),last_paths.end(),first_paths.begin(),first_paths.end(), 
      	  		      std::inserter(not_in_first, not_in_first.begin()));

      	  if(!not_in_last.empty() || !not_in_first.empty()) {

      	    //std::cout<< element.first << " : "  ;
      	    y -= pitch;  
	    y_x1.push_back(y);
	    s_x1.push_back(element.first);
    
      	    std::map<int,std::string> output; 
      	    int count(0);
      	    for (unsigned int iPath = 0; iPath < not_in_last.size(); ++iPath) {
      	      //std::cout << not_in_last[iPath] << " ; " ;
	      
	      if((output[count]+not_in_last[iPath]).length()<60){
		output[count]+=not_in_last[iPath];
		output[count]+=";";
	      } else {
		count++;
		output[count]+=not_in_last[iPath];
	      }
      	    }

      	    for (unsigned int br=0; br<output.size();br++){
	      y_x2.push_back(y-(br*pitch));
	      s_x2.push_back("#color[6]{"+output[br]+"}");
      	    }  
      	    //std::cout << " ||||||";
	    
      	    output.clear();
      	    int count1=0;
      	    for (unsigned int jPath = 0; jPath < not_in_first.size(); ++jPath) {
      	      //std::cout << not_in_first[jPath] << " ; " ;

	      if((output[count]+not_in_first[jPath]).length()<60){
		output[count]+=not_in_first[jPath];
		output[count]+=";";
	      } else {
		count++;
		output[count]+=not_in_first[jPath];
	      }
      	    }

      	    for (unsigned int br=0; br<output.size();br++){  
	      y_x3.push_back(y-(br*pitch));
	      s_x3.push_back("#color[8]{"+output[br]+"}");
      	    }

	    // decrease the y position to the maximum of the two lists
	    y-=std::max(count,count1)*pitch;
	    y_line.push_back(y-0.008);
	    
      	    //std::cout << std::endl;
 
      	  } // close if there is at least a difference 
      	} // if there is a common key
      }//loop on the keys

      TCanvas canvas("AlCaRecoTriggerBits","AlCaRecoTriggerBits",2500.,std::max(y_x1.size(),y_x2.size())*40); 

      TLatex l;
      // Draw the columns titles
      l.SetTextAlign(12);

      float newpitch = 1/(std::max(y_x1.size(),y_x2.size())*1.65);
      float factor  = newpitch/pitch;
      l.SetTextSize(newpitch-0.002);
      canvas.cd();
      for(unsigned int i=0;i<y_x1.size();i++){
	l.DrawLatexNDC(x1,1-(1-y_x1[i])*factor,s_x1[i].c_str());
      }

      for(unsigned int i=0;i<y_x2.size();i++){
	l.DrawLatexNDC(x2,1-(1-y_x2[i])*factor,s_x2[i].c_str());
      }

      for(unsigned int i=0;i<y_x3.size();i++){
	l.DrawLatexNDC(x3,1-(1-y_x3[i])*factor,s_x3[i].c_str());
      }

      canvas.cd();
      canvas.Update();

      TLine lines[y_line.size()];
      unsigned int iL=0;
      for (const auto & line : y_line){
	//std::cout<<1-(1-line)*factor<<std::endl;
	lines[iL] = TLine(gPad->GetUxmin(),1-(1-line)*factor,gPad->GetUxmax(),1-(1-line)*factor);
	lines[iL].SetLineWidth(1);
	lines[iL].SetLineStyle(9);
	lines[iL].SetLineColor(2);
	lines[iL].Draw("same");
	iL++;
      }

      //canvas.SetCanvasSize(2000,(1-y)*1000);
      std::string fileName(m_imageFileName);
      canvas.SaveAs(fileName.c_str());
      return true;
    }
  };

}


PAYLOAD_INSPECTOR_MODULE( AlCaRecoTriggerBits ){
  PAYLOAD_INSPECTOR_CLASS( AlCaRecoTriggerBits_Display );
  PAYLOAD_INSPECTOR_CLASS( AlCaRecoTriggerBits_Compare );
}
