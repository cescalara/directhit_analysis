#include <Riostream.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TColor.h>

//Stefano Durando

//The same of SHAPE_HidesAParticleEvent but it prints on a text file just the numbers of the GTUs of Direct Events, it prints this number in a column in order to let the analysis been easier.

void DirectEvents(const char* filein, const char* fileout, const int  maxGTU){
  
  
  ifstream in(filein);
  if(!in)
    {
      cout<<"The file "<<filein<<" does not exist "<<endl;
      return;
    }
  ofstream out(fileout);

  Int_t TH=6;

    TFile* f=new TFile(filein,"r");
    TTree* t= f -> Get("tevent");
    unsigned char pcd[1][1][48][48];
    t->SetBranchAddress("photon_count_data",pcd);
    bool  vectorevents[maxGTU];
    bool shape[maxGTU];
    Int_t eventscounts[maxGTU];

    // out<<endl<<"File: "<<filein<<endl<<endl<<"******************************"<<endl<<"threshold: "<<TH<<"counts"<<endl<<"******************************"<<endl;
    
    // out<<endl<<"Events:"<<endl<<endl;
    
    
    for(Int_t GTU=0;GTU<maxGTU;GTU++){
      
      t -> GetEntry(GTU);
      Int_t x;
      Int_t lstart=48,lend=0,hstart=48,hend=0;                  
      bool ve=0;                                          
      bool stop=0;
      Int_t* threshold=&TH;
      Int_t newTH=2*TH;
      Int_t count=0;
      
      for(Int_t row=0;row<48;row++){
	
	for(Int_t column=0;column<48;column++){
	  
	  x = pcd[0][0][column][row]; 
	  
	  if(x>2*TH)
	    threshold=&newTH;
	  
	  if(x>*threshold-*threshold/5){
	    
	    if(lstart>column)
	      lstart=column;
	    if(hstart>row)
	      hstart=row;
	    if(lend<column)
	      lend=column;
	    if(hend<row)
	      hend=row;
	  }
	  
	  if(x>*threshold){
	    
	    Int_t r=1,c=1,d=1;
	    if(x>count)
	      count=x;
	    
	    if(stop==0){
	      
	      while(r<48-row){
		
		Int_t u;
		u=pcd[0][0][column][row+r];
		
		if(u>*threshold-1-*threshold/5){
		  
		  r+=1;
		  
		  if(r>2){                                       
		    ve=1;
		    stop=1;
		  }
		}
		
		else
		  r=48-row;
		
	      }
	      
	      if(stop==0){
		
		Int_t s = 48-column;
		if(row>column)
		  s=48-row;
		
		while(d<s){
		  Int_t u=pcd[0][0][column+d][row+d];            
		  if(u>*threshold-1-*threshold/5){
		    d+=1;
		    if(d>2){
		      ve=1;
		      stop=1;
		    }
		  }
		  else
		    d=s;
		}
		
		if(stop==0){
		  
		  d=0;
		  Int_t ss =48-row;
		  if(ss>column)
		    ss=column;
		  
		  while(d<ss){
		    
		    Int_t u=pcd[0][0][column-d][row+d];
		    if(u>*threshold-1-*threshold/5){
		      d+=1;
		      if(d>2){
			ve=1;
			stop=1;
		      }
		    }
		    else
		      d=ss;
		  }
		  
		  
		  if(stop==0){
		    
		    while(c<48-column){
		      
		      Int_t u;
		      u=pcd[0][0][column+c][row];
		      if(u>*threshold-1-*threshold/5){
			c+=1;
			if(c>2){                                  
			  ve=1;
			}
		      }
		      else
			c=48-column;
		      
		    }     	 
		  }
		}
	      }
	    }
	  }
	}
      }
      
      threshold=&TH;
      
      if(ve==1){
	
	vectorevents[GTU]=1;
	Int_t  height=hend-hstart+1;
	Int_t  length=lend-lstart+1;

	eventscounts[GTU]=count;
	
	if(height/length<1/1.5||height/length>1.5)
	  shape[GTU]=1;
	if(length>30||height>30)
	  shape[GTU]=1;
	else
	  shape[GTU]=0;   
      }
      
      else
	vectorevents[GTU]=0;
    }
    
    for(Int_t GTU=0;GTU<maxGTU;GTU++){                 
      
      if(vectorevents[GTU]==1){
	Int_t i=0;
	while(vectorevents[GTU+i]==1)
	  i++;
	if(i>2){
	  while(i>=0){
	    vectorevents[GTU+i]=0;
	    i--;
	  }
	}
      }
    }
    
    float cluster=0,line=0;
    Int_t nbright=0;
    
    for(Int_t i=0; i<maxGTU; i++){
      if(vectorevents[i]==1)
	out<<i<<endl;    	  
    }	
    // 	out<<"{maximum count:"<<eventscounts[i]<<"}";
    // 	if(eventscounts[i]>20){
    // 	  out<<"SUPERBRIGHT";
    // 	  nbright+=1;
    // 	}
      // 	out<<endl;
      // }
    // }
   //  out<<endl;
//     float percentage=(cluster+line)/maxGTU;
//     float pline=line/(cluster+line);
//     float pcluster=cluster/(cluster+line);
    
//     out<<endl<<endl<<"*******************************************************************************************************************************"<<endl<<"# of GTUs: "<<maxGTU<<endl<<endl;
  
//     out<<"|# of direct events|   |% of direct events|   |% of LINES over the all DEs|   |% of CLUSTER over the all DEs|"<<endl<<endl;
//     out<<(cluster+line)<<"                       "<<percentage*100<<"%"<<"                       "<<pline*100<<"%"<<"                       "<<pcluster*100<<"%"<<endl;

//     out<<"% #Brightevents/#ofDirectevents: "<<endl<<"  "<<100*(nbright)/(cluster+line)<<"%"<<endl<<endl;
//     out<<endl<<"*******************************************************************************************************************************"<<endl;
//   }
}
