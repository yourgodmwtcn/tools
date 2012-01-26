/*
Read avisynth files into matlab( http://avisynth.org/mediawiki/Main_Page)
Image is returned as rgb

img = mavs( avisynthfile.avs, 1 );
infostruct = mavs( avisynthfile.avs );
mavs('cleanup')   releases the avisynth library & function pointers

An avisynth script can look as simple as this:
DirectShowSource( "avifile.avi" );

Warning, an avisynth file called e.g myfile.avs containing DirectShowSource("myfile.avs") seem to cause an infite loop


Compilation:
install Microsoft Visual C++ 2008
http://www.microsoft.com/Express/VC/#webInstall
"mex -setup" and specify this compiler
"mex mavs.cpp" to compile
*/


#include "windows.h"
#include "mex.h"
#include "avisynth\src\core\avisynth.h"
#include <iostream>
#include <string>
#include <vector>
#include <ctype.h>




// To compile:
// mex mavs.cpp


// define the type used for data samples 
#define TYPE UINT8_T 

//variables
static const unsigned int MEM_MAX = 64;		//maximum memory avisynth can use (MB)
static const unsigned int MAX_NAME_LENGTH = 1024; // maximum length of the filename (dll)
char FILENAME[MAX_NAME_LENGTH]; // Storage of the filename
static char AVISYNTH_DLL[] = "avisynth.dll";
static HINSTANCE AVISYNTH = NULL; // Handle to library (dll)
static IScriptEnvironment* (__stdcall *CreateEnv)(int) = NULL;
static IScriptEnvironment *env = NULL;
bool   info    = false; 
int    framenr = -1;
std::string ext;
static std::string currentFileName = "";
std::vector<std::string> args;
std::vector<std::string> currentArgs;
static PClip clip=0;
static const VideoInfo *vi;

//print debug info
const bool dbg = false;

//functions
static void usage();
static void cleanup();
void getClipAVSARGS( PClip &clip, std::vector<std::string> args );
void getClipAVS( PClip &clip, char * FILENAME );
void setRGB( PClip &clip, TYPE *imgOut, int h, int w );
PVideoFrame getFrame( PClip &clip );


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ){       
        
						
		//input    nrhs = number of inputs, nlhs is output
        if ( nrhs>2 || nrhs<1 ){
			usage();
		}	
		//output
		if ( nlhs>1 ){
			usage();
		}				
		
		int type = -1;
		args.clear();                        
		int n = mxGetNumberOfElements( prhs[0] );
		char * FILENAME = (char*) mxCalloc( n+1, sizeof(char) );			
        if ( mxIsChar( prhs[0]) ){
			//if first input is char, get filename					
			type = 0;
			// first input is filename
			if( n>=MAX_NAME_LENGTH ){
				mexErrMsgTxt( "Too large filename..." );
			}		
			mxGetString( prhs[0], FILENAME, n+1 );		
			std::string filename = FILENAME; 
			ext = filename.erase( 0, filename.size()-3 );			
		
		}else if( mxIsCell( prhs[0]) ){
			//if first input is cell, it`s a cell array of avisytnh arguments
			if(n==0){
				mexErrMsgTxt( "Empty cell array not allowed" );			
			}
			type = 1;
			//first input is a cell arary of avisynth arguments			
			for(int i=0;i<n;i++){
				std::string str = mxArrayToString( mxGetCell( prhs[0],i ) );
				args.push_back( str );
			}			
		}else{
			mexErrMsgTxt( "First input must be char(filename) or cell(avisynth args)" );
		}				
		
		info = false;
		if( nrhs==1 ){
			if( strcmp(FILENAME,"cleanup")==0 ){
				//if only one input, which equals "cleanup", cleanup...
				cleanup();
				return;
			}else{
				// if only one input (filename), return infostruct
				info = true;
			}
		}else if( nrhs==2 ){
			//second input is framenr
			if ( !mxIsDouble( prhs[1]) ){			
				mexErrMsgTxt( "Framenr (second argument) must be an positive integer" );			
			}
			framenr = (int)mxGetScalar( prhs[1] );
			if( framenr<0 ){
				mexErrMsgTxt( "Framenr (second argument) must be an positive integer" );
			}
			if(dbg){
				mexPrintf("Frame nr: %d\n",framenr);
			}
		}
		
		
        //------------------------------------------ Load library & functions(if not already loaded)
        if ( !AVISYNTH ) {
                AVISYNTH = LoadLibrary( TEXT(AVISYNTH_DLL) );
                if ( !AVISYNTH ) {
                        std::string errmsg = "Could not load library ";
                        errmsg += AVISYNTH_DLL;
                        mexErrMsgTxt( errmsg.c_str() );
                }
				if(dbg){
					mexPrintf("Avisynth library loaded\n");
					mexPrintf("Loading script environment CreateEnv\n");		
				}	
				CreateEnv = (IScriptEnvironment *(__stdcall *)(int)) GetProcAddress(AVISYNTH, "CreateScriptEnvironment");						
				if (!CreateEnv){
					std::string errmsg = "Could load avisynth function CreateScriptEnvironment ";
					mexErrMsgTxt( errmsg.c_str() );			
				}
				if(dbg){
					mexPrintf("Creating script environment env\n");		
				}
				env = CreateEnv(AVISYNTH_INTERFACE_VERSION);
				if (!env){
					std::string errmsg = "Could create avisynth script environment ";            
					mexErrMsgTxt( errmsg.c_str() );			
				}
				env->SetMemoryMax( MEM_MAX );
				mexAtExit( cleanup ); 
        }

        //------------------------------------------------finished loading   					
			bool clipChanged = false;
			switch(type){
			
				case 0:				
					// file is an avifile or avisynth script
					if( clip==0 || strcmp( currentFileName.c_str(), FILENAME )!=0 ){						
					// only set the clip object if the clip is a null pointer, or if the filename has changed (if using static clip)
						if(dbg){
							mexPrintf("FILENAME has changed\n");
						}
						std::string ext2;
						for(int i=0;i<(int)ext.size();i++){			
							ext2.push_back( (char)tolower(ext[i]) );
						}				
						if( strcmp( ext2.c_str(), "avs" )==0 ){
							if(dbg){
								mexPrintf("File is avs file\n");
							}
							//file is an avisynth script, try to evaluate it
							getClipAVS( clip, FILENAME );			
							clipChanged = true;
						}else{
							if(dbg){
								mexPrintf("File is movie file: %s\n",TEXT(FILENAME));
							}
							//file is a videofile, try to open it by using DirectShowSource("avifile.avi")
							args.push_back("DirectShowSource(\"");args[0] += TEXT(FILENAME);
							args[0] += "\")";												
							getClipAVSARGS( clip, args );
							clipChanged = true;
						}
					}
					currentFileName.assign( TEXT(FILENAME) );
					break;
				case 1:
					// try to evaluate cell array of avisynth arguments
					
					//check if arguments vector has changed
					bool equal; equal = false;	
					if( currentArgs.size()==args.size() && args.size()!=0 ){
						equal = true;
						for(int i=0;i<(int)args.size();i++){
							if( strcmp( args[i].c_str(), currentArgs[i].c_str() )!=0 ){
								equal = false;
								break;
							}
						}
					}
					if( !equal ){
						//only evaluate if arguments has changed
						if(dbg){
							mexPrintf("Arguments has changed\n");
						}
						currentArgs.clear();
						for(int i=0;i<(int)args.size();i++){
							currentArgs.push_back( args[i] );
						}						
						getClipAVSARGS( clip, args );
						clipChanged = true;
					}
					break;
				default:
					mexErrMsgTxt( "Could not determine type of input" );
			
		}
		
		if(dbg){		
			mexPrintf("Value of clip: %d\n",clip);
		}
		if(clip==0){
			mexErrMsgTxt( "Something went wrong..." );
		}
		
		if( clipChanged ){			
			vi = &clip->GetVideoInfo();
			if(dbg){
				mexPrintf("Height of clip: %d\n",vi->height);
			}
		}

		
		if( !info ){	
			//return frame
			
			if( framenr>vi->num_frames ){
				char errmsg[1024]; sprintf(errmsg, "Error, trying to access frame %d of %d", framenr, vi->num_frames );			
				mexErrMsgTxt( errmsg );
			}			
		
			const int h = vi->height; 
			const int w = vi->width; 			
			TYPE *imgOut = (TYPE*) mxCalloc( w*h*3, sizeof(TYPE) );						
			int dims[3]; dims[0]=h;dims[1]=w;dims[2]=3;
			plhs[0] = mxCreateNumericArray( 3, dims, mxUINT8_CLASS, mxREAL );			
			imgOut = (TYPE*)mxGetPr( plhs[0] );			
			setRGB( clip, imgOut, h ,w ); 
									
		}else{
			//return info struct			
			const char *infoFields[] = {"width","height","fps","numFrames"};
			plhs[0] = mxCreateStructMatrix( 1, 1, 4, infoFields );
			mxSetFieldByNumber(plhs[0],0,0,mxCreateDoubleScalar( vi->width ));
			mxSetFieldByNumber(plhs[0],0,1,mxCreateDoubleScalar( vi->height ));
			mxSetFieldByNumber(plhs[0],0,2,mxCreateDoubleScalar( vi->fps_numerator/vi->fps_denominator ));
			mxSetFieldByNumber(plhs[0],0,3,mxCreateDoubleScalar( vi->num_frames ));			
			
		}
		

}


static void cleanup(){
// Free library resources, and zero library/function pointers
        
		if(dbg){
			mexPrintf("Cleaning up\n");
		}
		if ( AVISYNTH ) {
			FreeLibrary(AVISYNTH);
        }
		AVISYNTH  = NULL;
		CreateEnv = NULL;
		env       = NULL;		
		//clip      = NULL;
		//clip = 0;
		
}

static void usage(){
        mexErrMsgTxt( "Usage:\n img  = mavs( 'FILENAME.avs', framenr )\n info = mavs( 'FILENAME.avs' )" );
}


PVideoFrame getFrame( PClip &clip ){

	if(dbg){
		mexPrintf("Getting frame\n");
	}
	try{
		PVideoFrame frame = clip->GetFrame( framenr, env );
		return frame;
	}catch( AvisynthError e ){
		mexErrMsgTxt( e.msg );			
		cleanup();		
	}		
	
}



void getClipAVSARGS( PClip &clip, std::vector<std::string> args ){
// try to open avifile using the avisynth arguments in args.

		if( args.size()==0 ){
			mexErrMsgTxt( "No avisynth arguments specified" );
		}		
		////mexPrintf("Trying to import video file\n");
		AVSValue res;
		for(int i=0;i<(int)args.size();i++){
			if(dbg){
				mexPrintf( "args[i] %s\n",args[i].c_str() );
			}
			try{
				res = env->Invoke( "Eval", args[i].c_str() );
				//res = env->Invoke( args[i].c_str(), AVSValue(&arg, 1));
			}catch( AvisynthError e ){
				mexErrMsgTxt( e.msg );
				cleanup();						
			}				
		}						
		if ( !res.IsClip() ){
			std::string errmsg = "Invalid avisynth file...";
            mexErrMsgTxt( errmsg.c_str() );		
		}
		//PClip *clip = new PClip();		
		clip = res.AsClip(); 
		//convert to rgb (nothing happens if already rgb)
		AVSValue args_conv[1] = { clip };
		clip = env->Invoke("ConvertToRGB24", AVSValue(args_conv, 1)).AsClip();												
		
}

void getClipAVS( PClip &clip, char * AVSFILE ){
//try to import avisynth file

		//////mexPrintf("Trying to import avisynth file\n");
		AVSValue args( TEXT(AVSFILE)  );		
		AVSValue res;
		try{
			res = env->Invoke("Import", AVSValue(&args, 1));
		}catch( AvisynthError e ){
			mexErrMsgTxt( e.msg );			
			cleanup();            
		}				
		if ( !res.IsClip() ){
			std::string errmsg = "Invalid avisynth file...";
            mexErrMsgTxt( errmsg.c_str() );		
		}
				
		if(dbg){
			mexPrintf("Trying to get clip\n");
		}
		clip = res.AsClip();
		
		//convert to rgb (nothing happens if already rgb)
		AVSValue args_conv[1] = { clip };
		clip = env->Invoke("ConvertToRGB24", AVSValue(args_conv, 1)).AsClip();			
		
}

void setRGB( PClip &clip, TYPE *imgOut, int h, int w ){
	// set the rgb values in the m*n*3 imgOut array
	
	if(dbg){
		mexPrintf("Setting rgb values\n");
	}
	PVideoFrame frame = clip->GetFrame( framenr, env );//getFrame( *clip );	
	
	//const int h = (frame)->GetHeight(); 		
	//const int w = (const int)((frame)->GetRowSize() / 3 );
	int pitch   = (frame)->GetPitch();	//pitch is the length of one line of data 	
	if(dbg){
		mexPrintf("h: %d  w: %d  pitch: %d\n",h,w,pitch);										
	}
		
	const unsigned char* srcp = (frame)->GetReadPtr(); //strange AviSynth RGB24 encoded data		
	int k=0;		
	for (int y=0; y<h; y++){ // Loop from bottom to top line (opposite of YUV colourspace).			
		k=0;	
		for (int x=0; x<pitch; x+=3) { // and from leftmost pixel to rightmost one.				
			//mexPrintf( "%d %d %d\n", (TYPE)srcp[x+2], (TYPE)srcp[x+1], (TYPE)srcp[x]);
			//mexPrintf(" %d\n", h - y+k*h);
			// h- is because rgb images are stored upside down in avisynth
			int ind1 = h - y+k*h -1;
			int ind2 = h - y+k*h + w*h -1;
			int ind3 = h - y+k*h + w*h*2 -1;			
			//if statements below are needed, because sometimes the pitch > w*3. Can lead to indexes out out range			
			if(ind1<w*h){
				imgOut[ ind1 ] = (TYPE)srcp[x+2];  //r
			}else{
				mexWarnMsgTxt ("Warning, file appears to be corrupt");	
			}
			if(ind2<w*h*2){
				imgOut[ ind2 ] = (TYPE)srcp[x+1];  //g
			}else{
				mexWarnMsgTxt ("Warning, file appears to be corrupt");	
			}
			if(ind3<w*h*3){
				imgOut[ ind3 ] = (TYPE)srcp[x];	   //b
			}else{
				mexWarnMsgTxt ("Warning, file appears to be corrupt");	
			}
			k++;
		}						
		srcp = srcp + pitch; // Add the pitch to the source pointer (go to next line).						
	}	
	
	
}

