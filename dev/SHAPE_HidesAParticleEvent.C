#include <Riostream.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TColor.h>

//Stefano Durando

//This macro Get as input the root file where data are stored and gives as output the numbers of Direct Events GTUs with Shape, maximum counts, size..

// filein is the root file where SPB data are stored.
// fileout is the text file where the output will be stored.
// maxGTU is the number of GTUs in the root file

void SHAPE_HidesAParticleEvent(const char* filein, const char* fileout, const int  maxGTU){
  
  
  ifstream in(filein);
  if(!in)
    {
      cout<<"The file "<<filein<<" does not exist "<<endl;
      return;
    }
  ofstream out(fileout);

  for(Int_t TH=6;TH<9;TH++){   //loop for different values of threshold

    TFile* f=new TFile(filein,"r");
    TTree* t= f -> Get("tevent");
    unsigned char pcd[1][1][48][48];
    t->SetBranchAddress("photon_count_data",pcd);
    bool  vectorevents[maxGTU];  // array through which storing the events (1 event, 0 none)
    bool shape[maxGTU];          // ... the shape between line and cluster.
    Int_t eventscounts[maxGTU]; // array for the maximum count of the matrix

    Int_t area[maxGTU];  // for the size ( with a really low accuracy)
    
    for(int i=0;i<maxGTU;i++)  
      area[i]=0;
    
    out<<endl<<"File: "<<filein<<endl<<endl<<"******************************"<<endl<<"threshold: "<<TH<<"counts"<<endl<<"******************************"<<endl;
    
    out<<endl<<"Events:"<<endl<<endl;

    
    // loop over the GTUs
    
    for(Int_t GTU=0;GTU<maxGTU;GTU++){ 
      
      t -> GetEntry(GTU);
      Int_t x;
      Int_t lstart=48,lend=0,hstart=48,hend=0;                  
      bool ve=0;                                          
      bool stop=0;
      Int_t* threshold=&TH;
      Int_t newTH=2*TH;
      Int_t count=0;

      //loop over each pixel of the matrix.
      
      for(Int_t row=0;row<48;row++){
	
	for(Int_t column=0;column<48;column++){
	  
	  x = pcd[0][0][column][row]; 
	  
	  if(x>2*TH)
	    threshold=&newTH; // if the event is too bright the threshold increases in order to define better the shape,..
	  
	  if(x>*threshold-*threshold/5){  //looks at the size of the events( height and length that will be calculated after).
	                                  
	    if(lstart>column)
	      lstart=column;
	    if(hstart>row)
	      hstart=row;
	    if(lend<column)
	      lend=column;
	    if(hend<row)
	      hend=row;
	  }
	  
	  // The event is considered catched (1 in vectorevents) if it has more than 2 adjacent pixels over threshold minus 20%Threshold in one of the four directions.
	  
	  if(x>*threshold){ //looks for events.
	    
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
		    ve=1; //event is catched if ve=1.
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
	area[GTU]=height*length;  //approximatively
	if(height/length<1/1.5||height/length>1.5)
	  shape[GTU]=1;
	if(length>30||height>30) // Direct Events larger than a quarter of the matrix are lines.
	  shape[GTU]=1;
	else
	  shape[GTU]=0;   
      }
      
      else
	vectorevents[GTU]=0;
    }

    //loops that select the event shorter than 3GTUs.
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
      if(vectorevents[i]==1){
	if(shape[i]==1){
	  out<<"GTU:"<<i<<" LINE , ";
	  line+=1;
	} 
	else{
	  out<<"GTU: "<<i<<" CLUSTER , ";
	  cluster+=1;
	}
	out<<"{maximum count:"<<eventscounts[i]<<"}, area: "<<area[i]<<" ";
	if(eventscounts[i]>20){ // to highlight brighter events.
	  out<<"SUPERBRIGHT";
	  nbright+=1;
	}
	out<<endl;
      }
    }
    out<<endl;
    float percentage=(cluster+line)/maxGTU;
    float pline=line/(cluster+line);
    float pcluster=cluster/(cluster+line);
    
    out<<endl<<endl<<"*******************************************************************************************************************************"<<endl<<"# of GTUs: "<<maxGTU<<endl<<endl;
  
    out<<"|# of direct events|   |% of direct events(over all the GTUs)|   |% of LINES over the all DEs|   |% of CLUSTER over the all DEs|"<<endl<<endl;
    out<<(cluster+line)<<"                       "<<percentage*100<<"%"<<"                       "<<pline*100<<"%"<<"                       "<<pcluster*100<<"%"<<endl;

    out<<endl<<"|Percentage of Brightevents/Directevents| "<<endl<<"      "<<100*(nbright)/(cluster+line)<<"%"<<endl<<endl;
    out<<endl<<"*******************************************************************************************************************************"<<endl;
  }
}
