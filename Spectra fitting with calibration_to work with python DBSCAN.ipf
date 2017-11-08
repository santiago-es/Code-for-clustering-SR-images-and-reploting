#pragma rtGlobals=3		// Use modern global access method and strict wave access.



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////// FUNCTION TO CALL MULTIPLE FILES //////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
20/05/2015

Coefficients:
              Estimate Std. Error t value Pr(>|t|)    
(Intercept) 16.8356709  0.9584892  17.565  1.1e-05 ***
clustx      -0.0017716  0.0009624  -1.841   0.1250    
clusty       0.0090973  0.0022610   4.024   0.0101 *  

2.31

// Main macro for loading experiments. 
Macro Multiple_Files(distance,firstframe,lastframe,xshift,crop,pixelwidth,calibpixel,intercept,xdrift,ydrift)
variable distance=240
Prompt distance, "Distance between localisation and centre of spectrum: "	
variable Firstframe=1
Prompt Firstframe, "First frame in the nile red section of the image: "
variable lastframe=4000
Prompt lastframe, "Last frame in the nile red section of the image: "				
variable xshift=-2
Prompt xshift, "Spectrum shift in X: "			
variable crop=350
Prompt crop, "Super-resolution image cropped by: "		
variable calibpixel=30
Prompt calibpixel, "Width used to fit spectra over in the calibrant: "	
variable pixelwidth=80
Prompt pixelwidth, "Width to fit spectra over: "
variable Intercept=16.8356709
Prompt intercept, "Intercept parameter from bead calibrant: "			
variable xdrift=-0.0017716
Prompt xdrift, "X-drift parameter from bead calibrant:"			
variable ydrift=0.0090973
Prompt ydrift, "Y-drift parameter from bead calibrant: "			
variable pixel=2.31
//Prompt pixel, "nM per pixel: "			

string imagename="Image.tif"								// Change the filenames here. 

string resultsname="Results3.csv"	
string clustername=""				// No longer needed

variable numberoffolders=5
make/o/n=(numberoffolders)/T filelist
make/o/n=1/T pathtosaveexperiment="BirdBox2:20170221_DRW_asyn_chaps_SPAINT:"
filelist[0]="BirdBox2:20170221_DRW_asyn_chaps_SPAINT:nr_fibrils:fibrils_1:"
filelist[1]="BirdBox2:20170221_DRW_asyn_chaps_SPAINT:nr_fibrils:fibrils_2:"
filelist[2]="BirdBox2:20170221_DRW_asyn_chaps_SPAINT:nr_fibrils:fibrils_3:"
filelist[3]="BirdBox2:20170221_DRW_asyn_chaps_SPAINT:nr_fibrils:fibrils_4:"
filelist[4]="BirdBox2:20170221_DRW_asyn_chaps_SPAINT:nr_fibrils:fibrils_5:"




Multirun(Imagename,Resultsname,Clustername,distance,xshift,crop,intercept,xdrift,ydrift,firstframe,lastframe,pixelwidth,calibpixel,pixel,numberoffolders)
end

	


function multirun(Imagename,Resultsname,Clustername,distance,xshift,crop,intercept,xdrift,ydrift,firstframe,lastframe,pixelwidth,calibpixel,pixel,numberoffolders)
string imagename,resultsname,clustername												// Various strings from the macro
variable distance,xshift,crop,intercept,xdrift,ydrift,firstframe,lastframe,pixelwidth,calibpixel,pixel,numberoffolders		// Various variables from the macro

wave/t pathtosaveexperiment
string saves=pathtosaveexperiment[0]
string saveto=saves+"Paramters.txt"
make/o/n=1 int=intercept,xd=xdrift,yd=ydrift
make/o/n=1 var_distance=distance,var_xshift=xshift,var_intercept=intercept,var_xdrift=xdrift,var_ydrift=ydrift,var_firstframe=firstframe,var_lastframe=lastframe,var_pixelwidth=pixelwidth,var_calibwidth=calibpixel,var_pixel=pixel		// Make the variables for later on

Save/J/W/O var_distance,var_xshift,var_intercept,var_xdrift,var_ydrift,var_firstframe,var_lastframe,var_pixelwidth,var_calibwidth,var_pixel as saveto				// Save the parameters to a text file for future reference. 



wave/t filelist
variable a
for(a=0;a<(numberoffolders);a+=1)				// Run through each of the folders. 
	
	string pathname=filelist[a]
	string folder=num2str(a)
	newmulti(folder,pathname,imagename,resultsname,clustername,distance,xshift,crop,intercept,xdrift,ydrift,pixel)
endfor
string saveto2=pathtosaveexperiment[0]+"Fitspectra.pxp"
saveexperiment as saveto2
end

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////RUN ON EACH FOLDER//////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

function newmulti(a,pathname,imagename,resultsname,clustername,distance,xshift,crop,intercept,xdrift,ydrift,pixel)
string a,pathname,imagename,resultsname,clustername
variable distance,xshift,crop,intercept,xdrift,ydrift,pixel
setdatafolder root:
newdatafolder/s $a
runmulti(pathname,imagename,resultsname,clustername,distance,xshift,crop,intercept,xdrift,ydrift,pixel)
end

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////LOAD AND FILTER OUT THE NON-CLUSTER FILES//////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

function runmulti(pathname,imagename,resultsname,clustername,distance,xshift,crop,intercept,xdrift,ydrift,pixel)
string pathname	,imagename,resultsname,clustername					// Various strings loaded
variable distance,xshift,crop,intercept,xdrift,ydrift,pixel		
string nametoload=pathname+imagename				
string fitstoload=pathname+"0:"+resultsname									
string clusttoload=pathname+clustername
print pathname
ImageLoad/T=tiff/S=0 /C=-1 /O /N=image nametoload					// Load the tiff image-original			
LoadWave/J/D/W/A/K=0 	fitstoload									// Load the fits
print fitstoload
wave cluster
print nametoload
make/o/n=(dimsize(cluster,0)) clusters
	variable ff
		for(ff=0;ff<(dimsize(clusters,0));ff+=1)
			clusters[ff]=cluster[ff]+1
		endfor
	
wave xw,yw				
duplicate/o xw,xpos														// Need to rename, since xw and yw are used elsewhere.
duplicate/o yw,ypos
//LoadWave/O/J/D/W/A/K=0 	clusttoload								// Load cluster details
wave xw
//duplicate/o xw,clusters													// Rename to cluster- ease of use. 
//wave clusters
variable num=(dimsize(xw,0))											// Total number of rows. 
wave Frame,OrigX,OrigY,Origvalue,error,noise,signal,snr,background,amplitude,angle,yw,x_sd,y_sd,precision,xpos,ypos
make/o/n=1 tFrame,tOrigX,tOrigY,tOrigvalue,terror,tnoise,tsignal,tsnr,tbackground,tamplitude,tangle,tyw,tx_sd,ty_sd,tprecision,txw,tclusters,txpos,typos		// These will contain the filtered values. 
variable a,b
for(a=0;a<num;a+=1)			// Run through each of the files. 
	if(clusters[a]==0)					
	else								// Only select those rows that are part of a cluster. 
	redimension/n=(b+1) tFrame,tOrigX,tOrigY,tOrigvalue,terror,tnoise,tsignal,tsnr,tbackground,tamplitude,tangle,tyw,tx_sd,ty_sd,tprecision,txw,tclusters,txpos,typos
	tframe[b]=frame[a]-1		// Remove 1 from each frame, since the imageJ indexing starts at 1, whereas igor indexing starts at 0. 
	torigx[b]=origx[a]
	torigy[b]=origy[a]
	torigvalue[b]=origvalue[a]
	terror[b]=error[a]
	tnoise[b]=noise[a]
	tsignal[b]=signal[a]
	tsnr[b]=snr[a]
	tbackground[b]=background[a]
	tamplitude[b]=amplitude[a]
	tangle[b]=angle[a]
	tyw[b]=yw[a]
	ty_sd[b]=y_sd[a]
	tx_sd[b]=x_sd[a]
	tprecision[b]=precision[a]
	tclusters[b]=clusters[a]
	txw[b]=xw[a]
	txpos[b]=xpos[a]
	typos[b]=ypos[a]
	b+=1
	endif
endfor
print b					
// Overwrite the original versions.
 
duplicate/o tframe,frame; KillWaves tframe
duplicate/o torigx,origx; KillWaves torigx
duplicate/o torigy,origy; KillWaves torigy
duplicate/o torigvalue,origvalue; KillWaves torigvalue
duplicate/o terror,error; KillWaves terror
duplicate/o tnoise,noise; KillWaves tnoise
duplicate/o tsignal,signal; KillWaves tsignal
duplicate/o tsnr,snr; KillWaves tsnr
duplicate/o tbackground,background; KillWaves tbackground
duplicate/o tamplitude,amplitude; KillWaves tamplitude
duplicate/o tangle,angle; KillWaves tangle
duplicate/o tyw,yw; KillWaves tyw
duplicate/o ty_sd,y_sd; KillWaves ty_sd
duplicate/o tx_sd,x_sd; KillWaves tx_sd
duplicate/o tprecision,precision; KillWaves tprecision
duplicate/o txw,xw; KillWaves txw
duplicate/o txpos,xpos; KillWaves txpos
duplicate/o typos,ypos; KillWaves typos


spectrum(distance,crop,intercept,xdrift,ydrift,xshift)
make/o/n=1 distance_t=distance,cropt_t=crop,intercept_t=intercept,xdrift_t=xdrift,ydrift_t=ydrift,xdrift_t=xdrift
fitspectra()
wave image
killwaves image
package(pathname)
accessandsave(pathname,pixel)
killdatafolder/z fits
end


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////This is to extract the spectra from the images//////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





function spectrum(distance,crop,intercept,xdrift,ydrift,xshift)		// Do fitting of Gaussians.
//Variables input from the call-to function
variable distance 											// Distance from localisation to the spectrum
variable crop 												// The original image is cropped- this needs to be inserted here.
variable intercept											// Drift intercept- from the calibration data from beads- essentially from the linear fit to ensure that positions and spectra match.
variable xdrift												// Drift x
variable ydrift												// Drift y
variable xshift												// If there is a drift in the x
wave superres,Source, Frame, End_Frame, origX, origY, origValue, Error, Noise, Signal, SNR, Background,ypos,xpos,Amplitude, Angle, XW, YW, X_SD, Y_SD, Precision		// Tell procedure that the waves loaded are needed.
String savedDF= GetDataFolder(1)
setdatafolder root:
wave var_pixelwidth,var_firstframe,var_lastframe,var_pixelwidth,var_calibwidth
variable pixelwidth=var_pixelwidth[0]
variable firstframe=var_firstframe[0]
variable calibpixel=var_calibwidth[0]
setdatafolder savedDF
//Following is to extract data from the loaded waves, and correct for in terms of the drift parameters, crop etc.
make/o/n=(dimsize(frame,0)) xcoord,ycoord,xcoordex,ycoordex, frames,signals,amplitudes		// Make waves to store the various data in
variable rows=(dimsize(frame,0))						// Get dimensions of the wave- useful for FOR loops. 
variable localisation										// Counter for FOR loop
for(localisation=0;localisation<rows;localisation+=1)		// Extract data from loaded wave.
	xcoord[localisation]=origx[localisation]
	ycoord[localisation]=origy[localisation]+crop
	xcoordex[localisation]=origx[localisation]			// Since cropped image used for localisations, then need to correct address of the localisation in the non-cropped image.
	ycoordex[localisation]=origy[localisation]+crop			// Since cropped image used for localisations, then need to correct address of the localisation in the non-cropped image.
	frames[localisation]=firstframe+frame[localisation]	
	signals[localisation]=signal[localisation]
	amplitudes[localisation]=amplitude[localisation]
	ypos[localisation]=ypos[localisation]+crop
endfor
// The fitting part is now to follow:
make/o/n=(dimsize(frame,0)) sums,centre,centrecorr,width,ampli,noise,rsq,intensity,sd,means,offset		// Make waves to store the various data in
wave image													// Need to call image.
variable i	

newdatafolder/o fits		
make/o/n=(pixelwidth) sumline								// Make a wave to store all of the data in.
variable halfwidth=round(pixelwidth/2)
for(i=0;i<(dimsize(frames,0));i+=1)
	make/o/n=(pixelwidth) templine											
	variable frametolookat=frames[i]							// Frame number	
	variable row=xcoord[i]+xshift										// column number. 
	variable colmin=ycoord[i]-distance-halfwidth
	variable colmax=ycoord[i]-distance+halfwidth
	variable j=0,k=0													// Populate line
		for(j=(colmin);j<(colmax);j+=1)
			templine[k]=((image[row-1][j][frametolookat]+image[row][row+1][frametolookat]+image[row+1][j][frametolookat])-(image[row-6][j][frametolookat]+image[row-5][row+5][frametolookat]+image[row+6][j][frametolookat]))/3		// Take average of three frames.
			sumline[k]=sumline[k]+templine[k]
			k+=1
		endfor
	
		string name=num2str(i)
		setdatafolder fits
		duplicate/o ::templine,$name		//Comment out to hide fits and save space. 	
		setdatafolder ::
endfor



end

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////Perform the fits//////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

function fitspectra()										// Function to fit spectra
wave sumline											// Sumline is all of the spectra combined
CurveFit/w=2/NTHR=0/TBOX=768 gauss  sumline/D 			// Fit the sum-line
wave w_coef
duplicate/o W_coef,guess								// Gueses for fits
wave xcoord											// Want to know number of points.
variable length=dimsize(xcoord,0)						// Determine length for number of points
setdatafolder fits										// Go into folder with the fits
duplicate/o ::guess,guess2								
variable a,b,c											// Running variables
make/o/n=(length) base,amp,xcent,width,chi,maxloc,maxi,lower,upper
for(a=0;a<(length);a+=1)								// Go through each of lines
string name=num2str(a)
duplicate/o $name,temp
K0 = guess2[0]/length;K1 = guess2[1]/length;K2 = guess2[2];K3 = guess2[3];				// Guesses- from fitting the sumline.
CurveFit/w=2/N=1/Q/G/NTHR=0/TBOX=768 gauss  temp /D 										// Fit the data.
wave W_coef										
base[a]=W_coef[0]
amp[a]=W_coef[1]																			// Store the data. 
xcent[a]=W_coef[2]					////////// CHANGED HERE			// If add to position, then subtract here. 
width[a]=W_coef[3]
chi[a]=V_chisq

variable n=100*(a/length)
print n
endfor
setdatafolder ::


end

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////See the particular spectrum//////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

function seespectra(num)					// See particular spectrum
variable num
wave sumline								// Sumline is all of the spectra combined
wavestats sumline							// Want to find maximum in sumline, since only want to fit the maximum in the spectrum, not the shoulders etc. 
variable centre_point=v_maxloc			// Location of the maximum

wave w_coef
duplicate/o W_coef,guess					// Gueses for fits
wave xcoord									// Want to know number of points.
variable length=dimsize(xcoord,0)			// Determine length for number of points
//variable length=10
setdatafolder fits
string graph=num2str(num)				// Display graph
duplicate/o $graph,temp
display temp
Label bottom "Pixel"
Label left "Intensity"
ModifyGraph mirror=1,minor=1,btLen=4,stLen=3
make/o/n=4 coeffs
wave base,amp,xcent,width,chi,lower,upper
coeffs[0]=base[num]
coeffs[1]=amp[num]
coeffs[2]=xcent[num]
coeffs[3]=width[num]
setdatafolder ::

k0=coeffs[0];k1=coeffs[1];k2=coeffs[2];k3=coeffs[3]
variable lowerbound=lower[num]
CurveFit/w=2/N=1/Q/G/NTHR=0/TBOX=768 gauss  temp /D 			// Without bounds
//CurveFit/G/H="1111"/NTHR=0 gauss  temp[lowerbound,upperbound] /D 

end

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////See the good spectrum//////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



function seespectragood(num2)					// See particular spectrum
variable num2
wave used
variable num=used[num2]
wave sumline								// Sumline is all of the spectra combined
wavestats sumline							// Want to find maximum in sumline, since only want to fit the maximum in the spectrum, not the shoulders etc. 
variable centre_point=v_maxloc			// Location of the maximum

wave w_coef
duplicate/o W_coef,guess					// Gueses for fits
wave xcoord									// Want to know number of points.
variable length=dimsize(xcoord,0)			// Determine length for number of points
//variable length=10
setdatafolder fits
string graph=num2str(num)				// Display graph
duplicate/o $graph,temp
display temp
Label bottom "Pixel"
Label left "Intensity"
ModifyGraph mirror=1,minor=1,btLen=4,stLen=3
make/o/n=4 coeffs
wave base,amp,xcent,width,chi,lower,upper
coeffs[0]=base[num]
coeffs[1]=amp[num]
coeffs[2]=xcent[num]
coeffs[3]=width[num]
setdatafolder ::

k0=coeffs[0];k1=coeffs[1];k2=coeffs[2];k3=coeffs[3]
variable lowerbound=lower[num]
CurveFit/w=2/N=1/Q/G/NTHR=0/TBOX=768 gauss  temp /D 			// Without bounds
//CurveFit/G/H="1111"/NTHR=0 gauss  temp[lowerbound,upperbound] /D 
end

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////See the bad spectrum//////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


function seespectrabad(num2)					// See particular spectrum
variable num2
wave notused
variable num=notused[num2]
wave sumline								// Sumline is all of the spectra combined
wavestats sumline							// Want to find maximum in sumline, since only want to fit the maximum in the spectrum, not the shoulders etc. 
variable centre_point=v_maxloc			// Location of the maximum

wave w_coef
duplicate/o W_coef,guess					// Gueses for fits
wave xcoord									// Want to know number of points.
variable length=dimsize(xcoord,0)			// Determine length for number of points
//variable length=10
setdatafolder fits
string graph=num2str(num)				// Display graph
duplicate/o $graph,temp
display temp
Label bottom "Pixel"
Label left "Intensity"
ModifyGraph mirror=1,minor=1,btLen=4,stLen=3
make/o/n=4 coeffs
wave base,amp,xcent,width,chi,lower,upper
coeffs[0]=base[num]
coeffs[1]=amp[num]
coeffs[2]=xcent[num]
coeffs[3]=width[num]
setdatafolder ::

k0=coeffs[0];k1=coeffs[1];k2=coeffs[2];k3=coeffs[3]
variable lowerbound=lower[num]
CurveFit/w=2/N=1/Q/G/NTHR=0/TBOX=768 gauss  temp /D 			// Without bounds
//CurveFit/G/H="1111"/NTHR=0 gauss  temp[lowerbound,upperbound] /D 
end

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////Combine all of the data//////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


function package(pathname)							// Make exportable for R cluster analysis
string pathname
setdatafolder fits
wave chi,width,xcent,amp,base,maxloc
duplicate/o chi,::m_chi
duplicate/o width,::m_width
duplicate/o xcent,::m_xcent
duplicate/o amp,::m_amp
duplicate/o base,::m_base
duplicate/o maxloc,::m_maxloc

setdatafolder ::
newdatafolder/o R_export

wave m_chi,m_width,m_xcent,m_amp,m_base,m_maxloc

make/o/n=((dimsize(xcent,0)),14) output
make/o/n=(1) used,notused
make/o/n=(dimsize(width,0)) use
variable a,b,c
wave frame,xcoordex,ycoordex,signals,centre,width,intensity,centrecorr,distance_t,cropt_t,intercept_t,xint_t,yint_t,xdrift_t,xcoord,ycoord,ydrift_t

for(a=0;a<(dimsize(m_chi,0));a+=1) // Go through all of data
if(sqrt(m_width[a]^2)<20 && sqrt(m_width[a]^2)>1.5)
if(m_xcent[a]>0)
if(m_amp[a]>0)
if(sqrt(centrecorr[a]^2)<30)
redimension/n=(b+1) used
use[a]=1
used[b]=a
b+=1
else
redimension/n=(c+1) notused
notused[c]=a
c+=1
endif
else
redimension/n=(c+1) notused
notused[c]=a
c+=1
endif
else
redimension/n=(c+1) notused
notused[c]=a
c+=1
endif
else
redimension/n=(c+1) notused
notused[c]=a
c+=1
endif
String savedDF= GetDataFolder(1)	// Remember CDF in a string.
setdatafolder root:
// Need these to convert to same pixel positions as the calibrant data. 
wave var_pixelwidth,var_calibwidth
variable pixelwidth=var_pixelwidth[0]
variable calibwidth=var_calibwidth[0]
variable halfpixel=round(pixelwidth/2	)			// This would usually equal 20
variable halfcalib=round(calibwidth/2	)			// This would usually equal 15
variable offset=halfpixel-halfcalib
SetDataFolder savedDF
output[a][0]=frame[a]
output[a][1]=xcoordex[a]
output[a][2]=ycoordex[a]
output[a][3]=signals[a]
output[a][4]=m_xcent[a]-offset														// Need to convert to centre to match the calibrant
output[a][5]=-offset+m_xcent[a]-(intercept_t[0]+xcoordex[a]*xdrift_t[0]+ycoordex[a]*ydrift_t[0])		// Correct for drift etc. and offset difference between calibrant and other.
output[a][6]=m_width[a]
output[a][7]=m_amp[a]
output[a][8]=m_base[a]
output[a][9]=m_amp[a]
output[a][10]=m_chi[a]


endfor
print b
print offset

setdatafolder R_export
duplicate/o ::output,output
duplicate/o ::use,use2
Movewave ::ypos,:yposi
Movewave ::Xpos,:xposi
Movewave ::xw,:xW
//duplicate/o ::output,output2

make/o/n=((dimsize(output,0))) frame,xcoordex,ycoordex,signal,centre,centrecorr,width,noise,ampli,chi,xpos,ypos
wave use
for(a=0;a<(dimsize(xpos,0));a+=1)
frame[a]=output[a][0]
xcoordex[a]=output[a][1]
ycoordex[a]=output[a][2]
signal[a]=output[a][3]
centre[a]=output[a][4]
centrecorr[a]=output[a][5]
width[a]=output[a][6]
ampli[a]=output[a][7]
noise[a]=output[a][8]
chi[a]=output[a][10]
endfor
killwaves output
wave path
setdatafolder ::
wave xW
string saveto=pathname+"FitSpectra.txt"
print saveto
setdatafolder R_export
wave chi,ampli,noise,width,centrecorr,centre,signal,ycoordex,xcoordex,frame,xW,use2,xposi,yposi
//Save/J/W/O chi,ampli,noise,width,centrecorr,centre,signal,ycoordex,xcoordex,frame,xW,use as saveto
Save/J/W/O chi,ampli,noise,width,centrecorr,centre,signal,ycoordex,xcoordex,xposi,yposi,frame,xW,use2 as saveto

setdatafolder ::
end


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////Calibrate functions below//////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

Macro Calibrate(distance,firstframe,lastframe,xshift,crop,pixelwidth)
variable distance=230
Prompt distance, "Distance between localisation and centre of spectrum: "	
variable Firstframe=1
Prompt Firstframe, "First frame in the nile red section of the image: "
variable lastframe=100
Prompt lastframe, "Last frame in the nile red section of the image: "				
variable xshift=0
Prompt xshift, "Spectrum shift in X: "			
variable crop=350
Prompt crop, "Super-resolution image cropped by: "			
variable pixelwidth=30
Prompt pixelwidth, "Width to fit spectra over: "

calibratefunction(distance,firstframe,lastframe,xshift,crop,pixelwidth)


end




function calibratefunction(distance,firstframe,lastframe,xshift,crop,pixelwidth)
variable distance,firstframe,lastframe,xshift,crop,pixelwidth	
make/o/n=1 var_distance=distance,var_firstframe=firstframe,var_pixelwidth=pixelwidth,var_lastframe=lastframe,var_xshift=xshift,var_crop=crop,var_calibwidth=pixelwidth
ImageLoad/T=tiff/S=0 /C=-1 /O /N=image 								// Load the tiff image-original
LoadWave/J/D/W/A/K=0 												// Load the fits
string path=s_path
wave xw,yw				
duplicate/o xw,xpos														// Need to rename, since xw and yw are used elsewhere.
duplicate/o yw,ypos
string pathname=path+":"
calibspectrum(distance,crop,xshift)
make/o/n=1 distance_t=distance,cropt_t=crop
fitspectra()
package(path)
wave image
killwaves image
string saveto2=path+"Fitspectra_long.pxp"
saveexperiment as saveto2
string saveto3=path+"Parameters_long.txt"
Save/J/W/O var_distance,var_firstframe,var_lastframe,var_xshift,var_crop,var_pixelwidth as saveto3

end

function calibspectrum(distance,crop,xshift)		// Do fitting of Gaussians.
//Variables input from the call-to function
variable distance 											// Distance from localisation to the spectrum
variable crop 												// The original image is cropped- this needs to be inserted here.
variable xshift												// If there is a drift in the x
wave superres,Source, Frame, End_Frame, origX, origY, origValue, Error, Noise, Signal, SNR, Background,ypos,xpos,Amplitude, Angle, XW, YW, X_SD, Y_SD, Precision		// Tell procedure that the waves loaded are needed.
String savedDF= GetDataFolder(1)
setdatafolder root:
wave var_pixelwidth,var_firstframe,var_lastframe,var_pixelwidth
variable pixelwidth=var_pixelwidth[0]
variable firstframe=var_firstframe[0]

setdatafolder savedDF
//Following is to extract data from the loaded waves, and correct for in terms of the drift parameters, crop etc.
make/o/n=(dimsize(frame,0)) xcoord,ycoord,xcoordex,ycoordex, frames,signals,amplitudes		// Make waves to store the various data in
variable rows=(dimsize(frame,0))						// Get dimensions of the wave- useful for FOR loops. 
variable localisation										// Counter for FOR loop
for(localisation=0;localisation<rows;localisation+=1)		// Extract data from loaded wave.
	xcoord[localisation]=origx[localisation]
	ycoord[localisation]=origy[localisation]+crop
	xcoordex[localisation]=origx[localisation]			// Since cropped image used for localisations, then need to correct address of the localisation in the non-cropped image.
	ycoordex[localisation]=origy[localisation]+crop			// Since cropped image used for localisations, then need to correct address of the localisation in the non-cropped image.
	frames[localisation]=firstframe+frame[localisation]	
	signals[localisation]=signal[localisation]
	amplitudes[localisation]=amplitude[localisation]
	ypos[localisation]=ypos[localisation]+crop
endfor
// The fitting part is now to follow:
make/o/n=(dimsize(frame,0)) sums,centre,centrecorr,width,ampli,noise,rsq,intensity,sd,means,offset		// Make waves to store the various data in
wave image													// Need to call image.
variable i	

newdatafolder/o fits		
make/o/n=(pixelwidth) sumline								// Make a wave to store all of the data in.
variable halfwidth=round(pixelwidth/2)
for(i=0;i<(dimsize(frame,0));i+=1)
	make/o/n=(pixelwidth) templine											
	variable frametolookat=frames[i]							// Frame number	
	variable row=xcoord[i]+xshift										// column number. 
	variable colmin=ycoord[i]-distance-halfwidth
	variable colmax=ycoord[i]-distance+halfwidth
	variable j=0,k=0													// Populate line
		for(j=(colmin);j<(colmax);j+=1)
			templine[k]=(image[row-1][j][frametolookat]+image[row][row+1][frametolookat]+image[row+1][j][frametolookat])/3		// Take average of three frames.
			sumline[k]=sumline[k]+templine[k]
			k+=1
		endfor
	
		string name=num2str(i)
		setdatafolder fits
		duplicate/o ::templine,$name		//Comment out to hide fits and save space. 	
		setdatafolder ::
endfor



end


macro see(spectrum)
variable spectrum=1
Prompt spectrum, "Spectrum to see "			// Set prompt for X param
seespectra(spectrum)
end

macro seegood(spectrum)
variable spectrum=1
Prompt spectrum, "Spectrum to see "			// Set prompt for X param
seespectragood(spectrum)
end


macro seebad(spectrum)
variable spectrum=1
Prompt spectrum, "Spectrum to see "			// Set prompt for X param
seespectrabad(spectrum)
end


function accessandsave(pathname,pixel)
string pathname
variable pixel



duplicate/o :R_export:xposi,xposi
duplicate/o :R_export:yposi,yposi
duplicate/o :R_export:centrecorr,centre_correct
//duplicate/o :R_export:yposi,yposi


wave xposi,yposi,x_sd,y_sd,background,noise,signal,frame,use,centre_correct,clusters

// Need to convert x-centers to wavelengths

variable a
duplicate/o centre_correct,wavelength

for(a=0;a<(dimsize(centre_correct,0));a+=1)
	wavelength[a]=581.5-centre_correct[a]*pixel
endfor

string saveto=pathname+"for_render.txt"
print saveto
Save/J/W/O frame,xposi,yposi,x_sd,y_sd,background,signal,wavelength,use as saveto

//for(c=0;c<(dimsize(xposi,0));c+=1)
	//xposi[c]=((xposi[c])/109)
	//yposi[c]=((yposi[c]-350)/109)+350
//endfor
	setdatafolder root:

end










