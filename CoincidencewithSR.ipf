#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.
function multiple()
variable numberoffolders=4
make/o/n=(numberoffolders)/T filelist

//filelist[0]="BirdBox2:20170414_Cell_exps:AB_bodies:1:561:0:"
//filelist[1]="BirdBox2:20170414_Cell_exps:AB_bodies:2:561:0:"
//filelist[2]="BirdBox2:20170414_Cell_exps:AB_bodies:3:561:0:"
//filelist[3]="BirdBox2:20170414_Cell_exps:AB_bodies:4:561:0:"
//filelist[4]="BirdBox2:20170414_Cell_exps:AB_bodies:5:561:0:"
//filelist[5]="BirdBox2:20170414_Cell_exps:AB_axons:1:561:0:"
//filelist[6]="BirdBox2:20170414_Cell_exps:AB_axons:2:561:0:"
//filelist[7]="BirdBox2:20170414_Cell_exps:AB_axons:3:561:0:"
filelist[0]="BirdBox2:20170414_Cell_exps:Media_2017-04-13_18-15-33:AB_axon:"
filelist[1]="BirdBox2:20170414_Cell_exps:Media_2017-04-13_18-15-33:AB_axoncontrol:"
filelist[2]="BirdBox2:20170414_Cell_exps:Media_2017-04-13_18-15-33:AB_body:"
filelist[3]="BirdBox2:20170414_Cell_exps:Media_2017-04-13_18-15-33:AB_bodycontrol:"


// Store data
make/o/n=(numberoffolders) coincidentclu,noncoincidentclu,fractionclu

variable f
//for(f=0;f<1;f+=1)
for(f=0;f<numberoffolders;f+=1)
	setdatafolder root:
	string folder=num2str(f)
	//PRINT FOLDER
	newdatafolder/s $folder
	
string path=filelist[f]
	
	loadSR(path)
	Loadtht(path)
	averageimages()
	//removebackground()
	//threshold()
	spots()
	coincident()
	
	countpoints()
	savewithheader(path)
	//countpoints()
	//xyextract()
	//distances()
	//plotdl(path)
	//savedlplot(path)
	
	
	stats(f)
setdatafolder root:
endfor
variable first=0,last=f
post_analyse(first,last)

end

//////////LOAD Super Resolution////////////


function loadSR(path)
string path
string load1=path+"561:3:Results1.csv"
print load1
LoadWave/J/D/W/K=0/A load1			


end


////////////LOAD ThT File///////////////////

function loadtht(path)
string path
string load1=path+"641:3:641.tif"

// Load - will put up prompt to look for image:
ImageLoad/T=tiff/S=0/C=500/LR3D /O /N=image load1

// Function to extract ThT and AF647 images

end

///////////////////STACK IMAGES////////////////////////

function averageimages()
wave image

// Duplicate ThT part of the image- in this case it's from frame 151-300:
duplicate/O/R=(0,512)(0,512)(1,200) image,ThT

// The ThT part of the image needs to be averaged over all of the frames:
imagetransform averageimage ThT

// Clean up unwanted images to save disk space:
killwaves tht
//killwaves image		 
//wave m_aveimage
//duplicate/o m_aveImage,Tht

wave m_stdvimage
duplicate/o m_stdvImage,Tht

//killwaves m_aveImage
end


// Function to remove backgrounds:

function removebackground()
wave tht
variable a,b,c,d,e,f
// ThT channel:

duplicate/o ThT,ThTfiltered,ThTbackground,fit_ThTbackground	// Make various waves of same size

	for(a=0;a<(dimsize(ThT,0));a+=1)			// Go through all of pixels
		for(b=0;b<(dimsize(ThT,1));b+=1)
			make/o temp=NAN				// Make a tempory section of the image
				for(d=-2;d<2;d+=1)
					for(f=-2;f<2;f+=1)			// This will take minimum pixel value over this range around each spot
						redimension/n=(e+1) temp
						temp[e]=tht[a+d][b+f]	
						e+=1
					endfor
				endfor
			e=0
			wavestats/q temp				
			thtbackground[a][b]=V_min		// Minimum pixel vvalue
			//thtfiltered[a][b]=tht[a][b]-V_min		// Remove background.
		endfor
	endfor
CurveFit/X=1/NTHR=0 Gauss2D  ThTbackground /D=fit_ThTbackground 
wave fit_ThTbackground


	for(a=0;a<(dimsize(ThT,0));a+=1)			// Go through all of pixels
		for(b=0;b<(dimsize(ThT,1));b+=1)
			
			thtfiltered[a][b]=tht[a][b]-fit_thtbackground[a][b]		// Remove background.
		endfor
	endfor
 // Remove unnecessary files:
 killwaves temp
 //killwaves ThTbackground


end

function threshold()
wave ThTfiltered

wavestats ThTfiltered

variable thresh=200
print thresh

make/o/n=(512,512) binary_image=0

variable a,b,c

for(a=0;a<512;a+=1)
	for(b=0;b<512;b+=1)
		if(ThTfiltered[a][b]>thresh)
			binary_image[a][b]=1
		endif
	endfor
endfor

end


//////////// OLD TYPE /////////////

function spots()
wave ThTFiltered,tht
duplicate/o tht,thtfiltered
// Do ThT Detection:	
variable patchsize =5									// If non-Gaussian model, then need a size of box to delete from working-image.
//variable thresholdtht =200									// Manual threshold.

imagestats ThTFiltered	
make/o/n=1 patchsizes=patchsize			
variable thresholdtht = 600								// Auto-threshold.
make/o/n=1 thtthreshold = thresholdtht								// For coincidence.
variable i,j,k,l,m
duplicate/o ThTFiltered,thtworking
make/o/n=((dimsize(thtfiltered,0)),(dimsize(thtfiltered,1))) binary_image=0
Make/O ThTX,ThTY,ThTIntensity,ThTxwidth,ThTywidth,ThTamp,ThTRatio		// Make various waves for data

// Thresholding gubbings:
do
redimension/n=(i+1) ThTX,ThTY,ThTIntensity,ThTxwidth,ThTywidth,ThTamp,thtratio		// Increase length of waves for new data. 	
	ImageStats  /M=1 thtworking					// Statistics of image																							
	ThTX[i]=V_maxRowLoc					// Maximum intensity stored.														
	ThTY[i]=V_maxColLoc
	ThTIntensity[i]=V_max	
																				
	variable row=V_maxRowLoc	
	variable col=V_maxColLoc


	
	
// This part is for straight-forward delete. 
	for(j=0;j<=(2*patchsize);j+=1)																								
		for(k=0;k<=2*patchsize;k+=1)
		ThTworking[V_maxRowLoc-patchsize+j][V_maxColLoc-patchsize+k][]=0
		binary_image[V_maxRowLoc-patchsize+j][V_maxColLoc-patchsize+k][]=1
		endfor
	endfor
	
	if(ThTx[i]==ThTx[i-1])				// Control to prevent double counting of spots. 
	print "ERROR"
	for(j=0;j<=(patchsize);j+=1)																								
		for(k=0;k<=patchsize;k+=1)
	ThTworking[row-patchsize+j][col-patchsize+k]=0
	
endfor
	endfor
	endif
	
	
	i+=1	
	
while (V_max>thresholdtht)
NewImage/K=0 Thtfiltered
ModifyImage Thtfiltered ctab= {*,*,Grays,0}
SetDrawLayer ProgFront
SetDrawEnv linefgc= (0,65535,0),fillpat= 0,xcoord= top,ycoord= left, save
// Part to draw rectangles:

for(i=0;i<(Dimsize(thtx,0));i+=1)
	
		DrawRect ThTX[i]-patchsize,ThTY[i]-patchsize,ThTX[i]+patchsize,ThTY[i]+patchsize
		variable a,b
	
endfor

end

/////////////////////


////////Check threshold

function checkthresh()
wave ThT,binary_image

newimage ThT
newimage binary_image

end



/////////Now need to check whether the aggregates are ThT active. 

function coincident()
wave binary_image,cluster,xw,yw

// First need to count how many clusters

wavestats/q cluster
variable numberofclusters=v_max
make/o/n=(v_max) coincident_clusters=0	// This is where to store whether they are clustered or not. 
make/o/n=(dimsize(xw,0)) coincidentpoints
variable a,b,c

for(a=0;a<(numberofclusters);a+=1)			// Go through each of the clusters
	make/o/n=1 tempx,tempy
		c=0
		for(b=0;b<(dimsize(xw,0));b+=1)
			if(cluster[b]==a)
				redimension/n=(c+1) tempx,tempy
				tempx[c]=xw[b]
				tempy[c]=yw[b]
				c+=1
			endif
		endfor
// Temp x and Tempy now contain co-ordinates of the cluster
	variable d
	variable coinc=0
	for(d=0;d<(dimsize(tempx,0));d+=1)
		variable xco=tempx[d]
		variable yco=tempy[d]
		if(binary_image[xco][yco]>0)
			coincident_clusters[a]+=1
			coinc=1
		endif
	endfor
	
	for(b=0;b<(dimsize(xw,0));b+=1)
			if(cluster[b]==a)
				coincidentpoints[b]=coinc
			endif
		endfor

endfor

make/o/n=(512,512) coincidentSR=0,nonCoincidentSR=0
wave origx,origy
for(a=0;a<(dimsize(xw,0));a+=1)
	if(coincidentpoints[a]>0 && cluster[a]>-1)
		xco=origx[a]
		yco=origy[a]
		coincidentsr[xco][yco]=1
		variable count
		count+=1
	elseif(cluster[a]>-1)
		xco=origx[a]
		yco=origy[a]
		noncoincidentsr[xco][yco]=1
	endif
endfor
wave tht
newimage ThT
newimage coincidentSR
newimage noncoincidentSR





end


function savewithheader(path)
string path

variable a,b,c,d
variable correct=0
//Load the files
make/o/n=1 all_Frame, all_origX, all_origY, all_origValue, all_Error, all_Noise, all_SNR,all_Background,all_Signal, all_Angle, all_XW, all_YW, all_X_SD, all_Y_SD, all_Precision__nm_
wave Source, Frame, origX, origY, origValue, Error, Noise, SNR,Background, Signal, Angle, XW, YW, X_SD,Y_SD, Precision__nm_,cluster,coincidentpoints
make/o/n=1 all_Framenoncoinc, all_origXnoncoinc, all_origYnoncoinc, all_origValuenoncoinc, all_Errornoncoinc, all_Noisenoncoinc, all_SNRnoncoinc,all_Backgroundnoncoinc,all_Signalnoncoinc, all_Anglenoncoinc, all_XWnoncoinc, all_YWnoncoinc, all_X_SDnoncoinc, all_Y_SDnoncoinc, all_Precision__nm_noncoinc
	
	for(b=0;b<(dimsize(frame,0));b+=1)
		if(cluster[b]>-1 && coincidentpoints[b]>0)
	redimension/n=(c+1) all_Frame, all_origX, all_origY, all_origValue, all_Error, all_Noise, all_SNR,all_Background,all_Signal, all_Angle, all_XW, all_YW, all_X_SD, all_Y_SD, all_Precision__nm_
		all_frame[c]=frame[b]+correct
		all_origX[c]=origx[b]
		all_origY[c]=origy[b]
		all_origValue[c]=origvalue[b]
		all_Error[c]=error[b]
		all_Noise[c]=noise[b]
		all_SNR[c]=snr[b]
		all_Background[c]=background[b]
		all_Signal[c]=signal[b]
		all_Angle[c]=angle[b]
		all_XW[c]=xw[b]
		all_YW[c]=yw[b]
		all_X_SD[c]=x_SD[b]
		all_Y_SD[c]=Y_sd[b]
		all_Precision__nm_[c]=Precision__nm_[b]
	
	c+=1
	else
	redimension/n=(d+1) all_Framenoncoinc, all_origXnoncoinc, all_origYnoncoinc, all_origValuenoncoinc, all_Errornoncoinc, all_Noisenoncoinc, all_SNRnoncoinc,all_Backgroundnoncoinc,all_Signalnoncoinc, all_Anglenoncoinc, all_XWnoncoinc, all_YWnoncoinc, all_X_SDnoncoinc, all_Y_SDnoncoinc, all_Precision__nm_noncoinc
		all_framenoncoinc[d]=frame[b]+correct
		all_origXnoncoinc[d]=origx[b]
		all_origYnoncoinc[d]=origy[b]
		all_origValuenoncoinc[d]=origvalue[b]
		all_Errornoncoinc[d]=error[b]
		all_Noisenoncoinc[d]=noise[b]
		all_SNRnoncoinc[d]=snr[b]
		all_Backgroundnoncoinc[d]=background[b]
		all_Signalnoncoinc[d]=signal[b]
		all_Anglenoncoinc[d]=angle[b]
		all_XWnoncoinc[d]=xw[b]
		all_YWnoncoinc[d]=yw[b]
		all_X_SDnoncoinc[d]=x_SD[b]
		all_Y_SDnoncoinc[d]=Y_sd[b]
		all_Precision__nm_noncoinc[d]=Precision__nm_[b]


	
	
	
	
	
	
	d+=1
	endif
	endfor
	
	wavestats/q all_frame
	correct=v_max




string tosave=path+"allfitswithheader_coincident.txt"



//Save/J/M="\n"/O/W all_Frame, all_origX, all_origY, all_origValue, all_Error, all_Noise,all_Background,all_Signal, all_Angle, all_XW, all_YW, all_X_SD, all_Y_SD, all_Precision__nm_ as tosave
variable f1 
Open f1 as tosave


fprintf f1, "#Localisation Results File\r#FileVersion Text.D0.E0.V2\r\r#Name Image (LSE)\r"
fprintf f1, "#Source <gdsc.smlm.ij.IJImageSource><singleFrame>0</singleFrame><extraFrames>0</extraFrames><path>/Volumes/BIRDBOX/20161012_DNAPAINT_Tau_Fids_2/1/01.tiff</path></gdsc.smlm.ij.IJImageSource>\r"
fprintf f1, "#Bounds x0 y0 w512 h512\r#Calibration <gdsc.smlm.results.Calibration><nmPerPixel>105.0</nmPerPixel><gain>55.5</gain><exposureTime>25.0</exposureTime><readNoise>0.0</readNoise><bias>500.0</bias><emCCD>false</emCCD></gdsc.smlm.results.Calibration>\r"
fprintf f1, "#Configuration <gdsc.smlm.engine.FitEngineConfiguration><fitConfiguration><fitCriteria>LEAST_SQUARED_ERROR</fitCriteria><delta>1.0E-4</delta><initialAngle>0.0</initialAngle><initialSD0>2.0</initialSD0><initialSD1>2.0</initialSD1><computeDeviations>"
fprintf f1, "false</computeDeviations><fitSolver>LVM</fitSolver><minIterations>0</minIterations><maxIterations>20</maxIterations><significantDigits>5</significantDigits><fitFunction>CIRCULAR</fitFunction><flags>20</flags><backgroundFitting>true</backgroundFitting>"
fprintf f1, "<notSignalFitting>false</notSignalFitting><coordinateShift>4.0</coordinateShift><signalThreshold>1665.0</signalThreshold><signalStrength>30.0</signalStrength><minPhotons>30.0</minPhotons><precisionThreshold>625.0</precisionThreshold><precisionUsingBackground>"
fprintf f1, "false</precisionUsingBackground><nmPerPixel>105.0</nmPerPixel><gain>55.5</gain><emCCD>false</emCCD><modelCamera>false</modelCamera><noise>0.0</noise><widthFactor>2.0</widthFactor><fitValidation>true</fitValidation><lambda>10.0</lambda><computeResiduals>false</computeResiduals>"
fprintf f1, "<duplicateDistance>0.5</duplicateDistance><bias>500.0</bias><readNoise>0.0</readNoise><maxFunctionEvaluations>1000</maxFunctionEvaluations><searchMethod>POWELL</searchMethod><gradientLineMinimisation>false</gradientLineMinimisation><relativeThreshold>1.0E-6</relativeThreshold>"
fprintf f1, "<absoluteThreshold>1.0E-16</absoluteThreshold></fitConfiguration><search>3.0</search><border>1.0</border><fitting>3.0</fitting><failuresLimit>10</failuresLimit><includeNeighbours>true</includeNeighbours><neighbourHeightThreshold>0.3</neighbourHeightThreshold><residualsThreshold>"
fprintf f1, "1.0</residualsThreshold><noiseMethod>QUICK_RESIDUALS_LEAST_MEAN_OF_SQUARES</noiseMethod><dataFilterType>SINGLE</dataFilterType><smooth><double>0.5</double></smooth><dataFilter><gdsc.smlm.engine.DataFilter>MEAN</gdsc.smlm.engine.DataFilter></dataFilter>"
fprintf f1, "</gdsc.smlm.engine.FitEngineConfiguration>\r"
fprintf f1, "#Frame\torigX\torigY\torigValue\tError\tNoise\tBackground\tSignal\tAngle\tX\tY\tXSD\tYSD\tPrecision"
wfprintf f1, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\r" all_Frame, all_origX, all_origY, all_origValue, all_Error, all_Noise,all_Background,all_Signal, all_Angle, all_XW, all_YW, all_X_SD, all_Y_SD, all_Precision__nm_
Close f1

string tosave2="Macintosh HD:Users:Mathew:Desktop:indexlist.txt"
variable f2
open/A f2 as tosave2
fprintf f2, "path[]=\"%s\"\r" path
close f2


tosave=path+"allfitswithheader_noncoincident.txt"



//Save/J/M="\n"/O/W all_Frame, all_origX, all_origY, all_origValue, all_Error, all_Noise,all_Background,all_Signal, all_Angle, all_XW, all_YW, all_X_SD, all_Y_SD, all_Precision__nm_ as tosave
 
Open f1 as tosave


fprintf f1, "#Localisation Results File\r#FileVersion Text.D0.E0.V2\r\r#Name Image (LSE)\r"
fprintf f1, "#Source <gdsc.smlm.ij.IJImageSource><singleFrame>0</singleFrame><extraFrames>0</extraFrames><path>/Volumes/BIRDBOX/20161012_DNAPAINT_Tau_Fids_2/1/01.tiff</path></gdsc.smlm.ij.IJImageSource>\r"
fprintf f1, "#Bounds x0 y0 w512 h512\r#Calibration <gdsc.smlm.results.Calibration><nmPerPixel>105.0</nmPerPixel><gain>55.5</gain><exposureTime>25.0</exposureTime><readNoise>0.0</readNoise><bias>500.0</bias><emCCD>false</emCCD></gdsc.smlm.results.Calibration>\r"
fprintf f1, "#Configuration <gdsc.smlm.engine.FitEngineConfiguration><fitConfiguration><fitCriteria>LEAST_SQUARED_ERROR</fitCriteria><delta>1.0E-4</delta><initialAngle>0.0</initialAngle><initialSD0>2.0</initialSD0><initialSD1>2.0</initialSD1><computeDeviations>"
fprintf f1, "false</computeDeviations><fitSolver>LVM</fitSolver><minIterations>0</minIterations><maxIterations>20</maxIterations><significantDigits>5</significantDigits><fitFunction>CIRCULAR</fitFunction><flags>20</flags><backgroundFitting>true</backgroundFitting>"
fprintf f1, "<notSignalFitting>false</notSignalFitting><coordinateShift>4.0</coordinateShift><signalThreshold>1665.0</signalThreshold><signalStrength>30.0</signalStrength><minPhotons>30.0</minPhotons><precisionThreshold>625.0</precisionThreshold><precisionUsingBackground>"
fprintf f1, "false</precisionUsingBackground><nmPerPixel>105.0</nmPerPixel><gain>55.5</gain><emCCD>false</emCCD><modelCamera>false</modelCamera><noise>0.0</noise><widthFactor>2.0</widthFactor><fitValidation>true</fitValidation><lambda>10.0</lambda><computeResiduals>false</computeResiduals>"
fprintf f1, "<duplicateDistance>0.5</duplicateDistance><bias>500.0</bias><readNoise>0.0</readNoise><maxFunctionEvaluations>1000</maxFunctionEvaluations><searchMethod>POWELL</searchMethod><gradientLineMinimisation>false</gradientLineMinimisation><relativeThreshold>1.0E-6</relativeThreshold>"
fprintf f1, "<absoluteThreshold>1.0E-16</absoluteThreshold></fitConfiguration><search>3.0</search><border>1.0</border><fitting>3.0</fitting><failuresLimit>10</failuresLimit><includeNeighbours>true</includeNeighbours><neighbourHeightThreshold>0.3</neighbourHeightThreshold><residualsThreshold>"
fprintf f1, "1.0</residualsThreshold><noiseMethod>QUICK_RESIDUALS_LEAST_MEAN_OF_SQUARES</noiseMethod><dataFilterType>SINGLE</dataFilterType><smooth><double>0.5</double></smooth><dataFilter><gdsc.smlm.engine.DataFilter>MEAN</gdsc.smlm.engine.DataFilter></dataFilter>"
fprintf f1, "</gdsc.smlm.engine.FitEngineConfiguration>\r"
fprintf f1, "#Frame\torigX\torigY\torigValue\tError\tNoise\tBackground\tSignal\tAngle\tX\tY\tXSD\tYSD\tPrecision"
wfprintf f1, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\r" all_Framenoncoinc, all_origXnoncoinc, all_origYnoncoinc, all_origValuenoncoinc, all_Errornoncoinc, all_Noisenoncoinc,all_Backgroundnoncoinc,all_Signalnoncoinc, all_Anglenoncoinc, all_XWnoncoinc, all_YWnoncoinc, all_X_SDnoncoinc, all_Y_SDnoncoinc, all_Precision__nm_noncoinc

Close f1

tosave2="Macintosh HD:Users:Mathew:Desktop:indexlist.txt"

open/A f2 as tosave2
fprintf f2, "path[]=\"%s\"\r" path
close f2




end



function stats(f)
variable f
wave coincident_clusters

variable a
variable c,nc
for(a=0;a<(dimsize(coincident_clusters,0));a+=1)
	if(coincident_clusters[a]>0)
		c+=1
	else
	nc+=1
	endif
endfor

variable fract=(c/(c+nc))

setdatafolder root:
wave coincidentclu,noncoincidentclu,fractionclu
coincidentclu[f]=c
noncoincidentclu[f]=nc
fractionclu[f]=fract




end














function kill()

	variable                      winMask;
 
	variable                      i,n;
	variable                      all=0x1000+0x40+0x10+0x4+0x2+0x1;
	string                        theWins;
 
	winMask = !winMask ? all : winMask;
 
	theWins = winList("*",";","WIN:"+num2iStr(winMask & all));
	for(i=0,n=itemsInList(theWins,";") ; i<n ; i+=1)
		doWindow/K $stringFromList(i,theWins,";");
	endfor;
end

macro close_windows()
kill()
end







////////////// DL ANALYSIS ETC.//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////// (1) First code to run from multirun /////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

function countpoints()							// How many points per cluster?
wave origX,origY,cluster								
variable a,b,c,z									// Some variables for later operations
wavestats cluster								// Number of clusters
variable cluster_max=V_max
make/o/n=(cluster_max) points_per_cluster		// Table to store counts	- looks at the highest cluster number in the table.


for(a=0;a<(cluster_max);a+=1)					
	c=0												// Reset the counter.
	for(b=0;b<(dimsize(cluster,0));b+=1)				// Go through all of cluster numbers
		if(cluster[b] ==(a+1))								// If the cluster number is equal to the number- add the +1 since we're not interested in cluster = 0. 
		c+=1											// Add 1 to c for each value that is equal to that cluster
		endif
	endfor
	points_per_cluster[a]=c
endfor

wavestats points_per_cluster					// Output the average and standard deviations of the counts per cluster. 
variable average_points_per_cluster=V_avg
variable sdev_points_per_cluster=V_sdev
string avg=num2str(average_points_per_cluster)
string std=num2str(sdev_points_per_cluster)


// order the clusters by size

make/o/n=(cluster_max) ordered_points,ordered_cluster_number
duplicate/o points_per_cluster,temp
variable d=0
	do
		
		redimension/n=(d+1) ordered_points,ordered_cluster_number
		ordered_points[d]=V_max
		ordered_cluster_number[d]=V_maxloc
		variable pos=v_maxloc
		temp[pos]=0
		wavestats/Q temp
		d+=1
	while (v_max>0)				// as long as expression is TRUE


end

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////Extract X and Y co-ordinates	//////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Needed to extract the x and y positions for each cluster. 
function XYExtract()								// Plot all clusters with colours corresponding to size.

wave ordered_cluster_number,cluster,origx,origY,xw,yw
duplicate/o xw,xp
duplicate/o yw,yp	
variable a,b,c,d,f
wave cluster,precision__nm_
duplicate/o precision__nm_,precision 
wavestats cluster
variable length=V_max							

for(a=0;a<(dimsize(ordered_cluster_number,0));a+=1)	// Looks for the highest value in the cluster.
	make/o/n=1 tempx,tempy,tempprec					// Need somewhere to store the data
	c=0
	variable cluster_number=ordered_cluster_number[a]+1 			// It has to be +1, because the ordered cluster number starts from zero, yet the cluster address wave starts from 1 (zero is a non-cluster address).
	for(b=0;b<(dimsize(cluster,0));b+=1)				// Go through all of the localisations information
		if(cluster[b]==cluster_number)
			variable x=xp[b]								// Set X and Y variables- from xp and yp which are the exact co-ordinates. 
			variable y=yp[b]
			variable sig=precision[b]					// Precision
			redimension/n=(c+1) tempx,tempy,tempprec		// Make wave longer for the files. 
			tempx[c]=Xp[b]							// Store data into tempx and tempy
			tempy[c]=yp[b]
			tempprec[c]=sig
			c+=1
		endif
	endfor
	string xval="x"+num2str(a)
	string yval="y"+num2str(a)
	string precval="prex"+num2str(a)
	duplicate/o tempx,$xval
	duplicate/o tempy,$yval
	duplicate/o tempprec,$precval
endfor
killwaves tempx,tempy,tempprec
end


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////// DL Plots	////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///THIS FOLLOWS THE CLUSTER ANALYSIS DONE IN THE PAPER///

// that it is 40 nm for now (which is ~0.3 pixels).

function distances()
wave cluster,origx,origy,xp,yp			
variable a,b,c,d,e,f
make/o/n=1 alldistances
wavestats cluster							// This is on the clusters- calculate how many clusters there are. 
variable length=v_max						

variable sig=0.3
make/o/n=(length) Reff			// This is the effective resolution- this will be equal to sqrt(rnn^2+sigma^2), where rnn = mean n.n. distance, and sigma=mean prec
								// Gould et al 2009.
make/o/n=(length) DLdist		// The neighbour distance

// Molecular positions need to be colour coded in terms of local density, where the local density is defined as the number of molecules
// within 5* the mean nearest neighbour distanceof all molecules in that cluster.

// First step- find out the mean nearest neighbour distance for all molecules in each PSD:

for(a=0;a<length;a+=1) 					// Go through each of clusters
	string xvar="X"+num2str(a)			// These are the cluster waves that need using. 
	string yvar="Y"+num2str(a)
	string precvar="prex"+num2str(a)
	make/o/n=1 tempdist				// Temporary wave in which to store the distances. 
	duplicate/o $xvar,tempx				// Contains all of x-coords.
	duplicate/o $yvar,tempy				// Contains all of y-coords
	duplicate/o $precvar,tempprec		// Contains all of precisions
	make/o/n=(dimsize(tempy,0)) tempcartdistance,tempNN
	variable dd
		for(b=0;b<(dimsize(tempx,0));b+=1)
 			 variable cartdistance=10000000		// This is just a large start value- it will immediately be replaced by the first distance. 
				for(c=0;c<(dimsize(tempx,0));c+=1)
	  				if(b==c)				// Don't want to measure the distance between a point and itself- this would = 0
	 				else 
						variable old_cartdistance=cartdistance		// Set the variable of the old-cartdistance to the previous distance measured. 
						variable xlength=tempx[c]-tempx[b]			// Calculate the distances between the points. 
						variable ylength=tempy[c]-tempy[b]
						cartdistance=sqrt(xlength^2+ylength^2)		// Pythag to calculate cartesian distance. 
							if(cartdistance<old_cartdistance)		// If the newly calculated cartesian distance is less than the old one	
																	// do nothing. 
							else				
							cartdistance=old_cartdistance			// Else let the cartesian distance equal the old one. I.e. If the cart distance is longer, then replace it with
																	// the one from the previous run. By doing this we will get to the neirest neighbour distance eventually. 
							endif
					endif
				endfor
			tempcartdistance[b]=cartdistance						// This is the nearest neighbour distance for each point in that cluster now. 
		endfor
	string toname2="Min"+num2str(a)								// Just need to copy this to its own wave now. 
	duplicate/o tempcartdistance,$toname2
	wavestats/q tempcartdistance									// Perform statisitics on it to get the average. 
	variable thresh=5*v_avg										// This is the threshold distance 
	dldist[a]=thresh													// Store this distance in the dldistance wave.
	variable tempcartnm=109*v_avg
	wavestats/q tempprec											// Get stats on the s.r. precisions.
	
	variable median_prec=median(tempprec)
	variable median_cart=109*median(tempcartdistance)	
	//variable reffective=sqrt(tempcartnm^2+v_avg^2)				// As defined in paper- the effective resolution. 
	variable reffective=sqrt(median_prec^2+median_cart^2)
	reff[a]=reffective												// Store in the wave. 		Gould et al. 2009. 
	
	// We now need to calculate how many neighbours each point has within a range of the threshold distance. 
	for(b=0;b<(dimsize(tempx,0));b+=1)								// Go through each of the points. 
		for(c=0;c<(dimsize(tempx,0));c+=1)							// Go through each of the points again to compare with the previous. 
			if(b==c)													// Don't compare with self.
			else
				xlength=tempx[c]-tempx[b]							// Calculate the lengths. 
				ylength=tempy[c]-tempy[b]
				cartdistance=sqrt(xlength^2+ylength^2)				// Calculate the cartesian distance.
				if(cartdistance<thresh)
					dd+=1											// If within the radius, then add 1 to variable dd. 
				endif
				endif
		endfor
		tempNN[b]=dd								// Save into temporary file.
	dd=0		// Reset variable 
	endfor
	string toname="DL"+num2str(a)				
	duplicate/o tempNN,$toname					// Rename as required. This cotains information about all of the nearest neighbours. 
endfor

end


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////// Make pretty images	/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

function plotdl(path)
string path
// Need to make an SR image that Matt likes. 
variable scale = 8									// Set the scale
variable pixelwidth=109								// nm pixel
make/o/n=((512*scale),(512*scale)) DL_Image=0,L_image=0,NN_Image=0		// Make the image- the dimensions may need to be altered, depending on the size of the image (this should be 10x the dimensions of the image)
wave xw
wavestats xw
variable length=v_max
variable a,b,c
make/o/n=1 xcoords,ycoords,NumNN,PrecisionList,NNdist,Clusternumber
for(a=0;a<length;a+=1)
	string xlist="x"+num2str(a)			// xcoord
	string ylist="y"+num2str(a)			// ycoord
	string dlist="dl"+num2str(a)			// no of nns
	string preclist="prex"+num2str(a)	// precision for sr fit
	string minlist="min"+num2str(a)		// nearest neighbour distance
	
	duplicate/o $xlist,tempx
	duplicate/o $ylist,tempy
	duplicate/o $dlist,tempd
	duplicate/o $preclist,tempprec
	duplicate/o $minlist,tempmin
		for(b=0;b<(dimsize(tempx,0));b+=1)
			redimension/n=(c+1) xcoords,ycoords,NNdist,PrecisionList,NumNN,Clusternumber
			
			xcoords[c]=tempx[b]
			ycoords[c]=tempy[b]
			NumNN[c]=tempd[b]
			PrecisionList[c]=tempprec[b]
			NNdist[c]=tempmin[b]*109
			clusternumber[c]=a+1
			c+=1
		endfor
endfor
string tosave=path+"All_Localisation_Information.txt"
Save/O/J/M="\n"/W Clusternumber,xcoords,ycoords,NNdist,PrecisionList,NumNN  as tosave


variable aa

for(a=0;a<(dimsize(xcoords,0));a+=1)							// Populate the matrix
	variable xpos=round(scale*(xcoords[a]))							// Get the co-ordinates
	variable ypos=round(scale*(ycoords[a]))		
					
	DL_image[xpos][ypos]=DL_image[xpos][ypos]+NumNN[a]							// Populate the matrix
	L_image[xpos][ypos]+=1
endfor



for(a=0;a<(dimsize(dl_image,0));a+=1)
	for(b=0;b<(dimsize(dl_image,1));b+=1)
		if(L_image[a][b]>0)
			NN_image[a][b]=DL_image[a][b]/L_image[a][b]
		else
			NN_image[a][b]=0
		endif
	endfor
endfor

wavestats/q xcoords
variable xmin=v_min
variable xmax=v_max

wavestats/q ycoords
variable ymin=v_min
variable ymax=V_max


wave dldist,reff,points_per_cluster

duplicate/o DLdist,DL_Distance
for(a=0;a<(dimsize(dldist,0));a+=1)
DL_Distance[a]=109*DL_Distance[a]
endfor
string tosave2=path+"All_Cluster_Information_all.txt"
Save/O/J/M="\n"/W points_per_cluster,Reff,DL_Distance  as tosave2



string image1=path+"NN_Image_all.tif"
string image2=path+"SR_Image_all.tif"

ImageSave/O/F/T="TIFF" NN_Image as image1
ImageSave/O/F/T="TIFF" L_Image as image2


end



macro plot()
kill()
plotall()
clusterplotpaper()
end

function plotall()
variable clustnumber
wave NN_image
NewImage/K=0 NN_Image
ModifyImage NN_Image ctab={*,*,Rainbow,0}
ModifyImage NN_Image ctab= {*,*,Rainbow,1}
SetDrawLayer ProgFront
SetDrawEnv linefgc= (65535,65535,65535),fillpat= 0,xcoord= top,ycoord= left, save

variable c
wave xw
wavestats xw
variable maxi=v_max
for(c=0;c<(maxi);c+=1)
string xlist="x"+num2str(c)
string ylist="y"+num2str(c)
string dllist="dl"+num2str(c)
duplicate/o $xlist,tempx
duplicate/o $ylist,tempy
wavestats/q tempx
//variable xmin=round(v_min-100)*8
//variable xmax=round(v_max-100)*8
//variable xdiff=round(v_max-v_min)

variable xmin=round(v_min)*8
variable xmax=round(v_max)*8
variable xdiff=round(v_max-v_min)

wavestats/q tempy
//variable ymin=round(v_min-100)*8
//variable ymax=round(v_max-100)*8
//variable ydiff=round(v_max-v_min)

variable ymin=round(v_min)*8
variable ymax=round(v_max)*8
variable ydiff=round(v_max-v_min)

DrawRect xmin,ymin,xmax,ymax




string clu=num2str(c+1)
SetDrawEnv textrgb= (65535,65535,65535),fsize= 10;DelayUpdate
DrawText xmax,ymax,clu
//print xmin
//print xmax
endfor


end


function clusterplot()
variable clustnumber
	Prompt clustnumber, "Cluster to plot: "		// Set prompt for x param
	
	DoPrompt "Cluster", clustnumber
	if (V_Flag)
		return -1								// User canceled
	endif

	
string xlist="x"+num2str(clustnumber)
string ylist="y"+num2str(clustnumber)
string dllist="dl"+num2str(clustnumber)
duplicate/o $xlist,tempx
duplicate/o $ylist,tempy
duplicate/o $dllist,tempd
variable scale=8
wavestats tempx
variable xmin=round(v_min)
variable xmax=round(v_max)
variable xdiff=round(v_max-v_min)

wavestats tempy
variable ymin=round(v_min)
variable ymax=round(v_max)
variable ydiff=round(v_max-v_min)

make/o/n=(scale*(xdiff+2),scale*(ydiff+2)) clu_DLimage=0,clu_Limage=0,clu_NNimage=0
variable a
for(a=0;a<(dimsize(tempx,0));a+=1)							// Populate the matrix
	variable xpos=round(scale*(tempx[a]-xmin+1))							// Get the co-ordinates
	variable ypos=round(scale*(tempy[a]-ymin+1))		
					
	clu_DLimage[xpos][ypos]=Clu_DLimage[xpos][ypos]+tempd[a]							// Populate the matrix
	clu_Limage[xpos][ypos]+=1
endfor

variable b

for(a=0;a<(dimsize(clu_dlimage,0));a+=1)
	for(b=0;b<(dimsize(clu_dlimage,1));b+=1)
		if(Clu_Limage[a][b]>0)
			Clu_NNimage[a][b]=Clu_DLimage[a][b]/clu_Limage[a][b]
		else
			clu_NNimage[a][b]=0
		endif
	endfor
endfor

Display;AppendMatrixContour clu_NNimage
ModifyContour clu_NNimage labels=0,autoLevels={*,*,100}
ModifyGraph mode=3,marker=19,msize=1
ModifyContour clu_NNimage autoLevels={*,*,10}
ModifyContour clu_NNimage autoLevels={*,*,50}
ModifyGraph mode=3,marker=19,msize=2
ColorScale/C/N=text0/A=MC  ctab={0,100,Rainbow,0}
ColorScale/C/N=text0/E
ColorScale/C/N=text0/A=MB/X=0.00/Y=0.00 vert=0
ModifyContour clu_NNimage ctabLines={*,*,Rainbow,1}
newimage clu_NNimage
ModifyImage clu_NNimage ctab={*,*,Rainbow,0}
ModifyImage clu_NNimage ctab= {*,*,Grays,0}
ImageThreshold/T=40 clu_nnimage
wave M_imagethresh
newimage M_imagethresh
ImageAnalyzeParticles /E/W/Q/M=3/A=10 stats,M_imageThresh
  end
  macro plotallnice()
  
  combineall()
  clusterplotpaperall()
  end
  
function combineall()
make/o/n=1 allx,ally,alldl

variable a,b,c,d
wave xw

wavestats xw
variable tot=v_max

for(a=0;a<tot;a+=1)

string xlist="x"+num2str(a)
string ylist="y"+num2str(a)
string dllist="dl"+num2str(a)

duplicate/o $xlist,tempx
duplicate/o $ylist,tempy
duplicate/o $dllist,tempd
	for(b=0;b<(dimsize(tempx,0));b+=1)
		redimension/n=(c+1) allx,ally,alldl
		
		allx[c]=tempx[b]
		ally[c]=tempy[b]
		alldl[c]=tempd[b]
		c+=1
	endfor
endfor

end
 
  
  
  
 function clusterplotpaperall()
 
variable upper=200

	Prompt Upper, "Upper value for scale: "		// Set prompt for x param
	
	DoPrompt "Upper", upper
	if (V_Flag)
		return -1								// User canceled
	endif

wave allx,ally,alldl
duplicate/o allx,tempx
duplicate/o ally,tempy
duplicate/o alldl,tempd

display tempx vs tempy

// Need to go through and plot a different color for each different density- i.e. 1 NN - Purple 100 NN - Red

// Get colors wave:
ColorTab2Wave Rainbow							// Make colour table. 

wave m_colors
variable increment=upper/100 
display
variable a,b,e=0
for(a=0;a<upper;a+=increment)
	string xcol="xco"+num2str(a)
	string ycol="yco"+num2str(a)
	make/o/n=1 tempxc,tempyc
	variable c=0
		for(b=0;b<(dimsize(tempd,0));b+=1)
			if(tempd[b]<(a+increment)&&tempd[b]>(a))
				redimension/n=(c+1) tempxc,tempyc
					tempxc[c]=tempx[b]
					tempyc[c]=tempy[b]
					c+=1
			endif
		endfor
		
		
		
		
	duplicate/o tempxc,$xcol
	duplicate/o tempyc,$ycol
	appendtograph $ycol vs $xcol
	
	
	variable r=m_colors[99-e][0]
	variable g=m_colors[99-e][1]
	variable bl=m_colors[99-e][2]
	ModifyGraph rgb($ycol)=(r,g,bl)
	
	e+=1
endfor
ModifyGraph mode=3,marker=19,msize=2
ModifyGraph msize(yco19)=3
ModifyGraph width=566.929,height=566.929
ModifyGraph mirror=1,minor=1,btLen=4,stLen=3
ColorScale/C/N=text1/F=0/A=RC/E  ctab={0,upper,Rainbow,1}
end


function clusterplotpaper()
variable clustnumber
	Prompt clustnumber, "Cluster to plot: "		// Set prompt for x param
	
	DoPrompt "Cluster", clustnumber
	if (V_Flag)
		return -1								// User canceled
	endif

	
string xlist="x"+num2str(clustnumber)
string ylist="y"+num2str(clustnumber)
string dllist="dl"+num2str(clustnumber)
duplicate/o $xlist,tempx
duplicate/o $ylist,tempy
duplicate/o $dllist,tempd

display tempx vs tempy

// Need to go through and plot a different color for each different density- i.e. 1 NN - Purple 100 NN - Red

// Get colors wave:
ColorTab2Wave Rainbow							// Make colour table. 

wave m_colors

display
variable a,b,e=0
for(a=0;a<100;a+=1)
	string xcol="xco"+num2str(a)
	string ycol="yco"+num2str(a)
	make/o/n=1 tempxc,tempyc
	variable c=0
		for(b=0;b<(dimsize(tempd,0));b+=1)
			if(tempd[b]<(a+2)&&tempd[b]>(a))
				redimension/n=(c+1) tempxc,tempyc
					tempxc[c]=tempx[b]
					tempyc[c]=tempy[b]
					c+=1
			endif
		endfor
		
		
		
		
	duplicate/o tempxc,$xcol
	duplicate/o tempyc,$ycol
	appendtograph $ycol vs $xcol
	
	
	variable r=m_colors[99-e][0]
	variable g=m_colors[99-e][1]
	variable bl=m_colors[99-e][2]
	ModifyGraph rgb($ycol)=(r,g,bl)
	e+=1
endfor
ModifyGraph mode=3,marker=19,msize=2
ModifyGraph msize(yco19)=3
end




  
 
function Particle_analysis(path)
string path
make/o/n=1 major,minor,particleclusternumber
wave xw
variable d
wavestats xw
variable length=v_max
make/o/n=(length) numberofparticles
variable l
for(l=0;l<length;l+=1)
variable clustnumber=l
kill()
string xlist="x"+num2str(clustnumber)
string ylist="y"+num2str(clustnumber)
string dllist="dl"+num2str(clustnumber)
duplicate/o $xlist,tempx
duplicate/o $ylist,tempy
duplicate/o $dllist,tempd
variable scale=8
wavestats/q tempx
variable xmin=round(v_min)
variable xmax=round(v_max)
variable xdiff=round(v_max-v_min)

wavestats/q tempy
variable ymin=round(v_min)
variable ymax=round(v_max)
variable ydiff=round(v_max-v_min)

make/o/n=(scale*(xdiff+2),scale*(ydiff+2)) clu_DLimage=0,clu_Limage=0,clu_NNimage=0
variable a
for(a=0;a<(dimsize(tempx,0));a+=1)							// Populate the matrix
	variable xpos=round(scale*(tempx[a]-xmin+1))							// Get the co-ordinates
	variable ypos=round(scale*(tempy[a]-ymin+1))		
					
	clu_DLimage[xpos][ypos]=Clu_DLimage[xpos][ypos]+tempd[a]							// Populate the matrix
	clu_Limage[xpos][ypos]+=1
endfor

variable b

for(a=0;a<(dimsize(clu_dlimage,0));a+=1)
	for(b=0;b<(dimsize(clu_dlimage,1));b+=1)
		if(Clu_Limage[a][b]>0)
			Clu_NNimage[a][b]=Clu_DLimage[a][b]/clu_Limage[a][b]
		else
			clu_NNimage[a][b]=0
		endif
	endfor
endfor

ImageThreshold/T=30/I clu_nnimage

wave M_imagethresh

ImageAnalyzeParticles/E/W/Q/M=1/A=5 stats,M_imageThresh
//newimage clu_nnimage
wave M_moments
wave w_boundaryy,w_boundaryx


for(a=0;a<dimsize(m_moments,0);a+=1)
redimension/n=(d+1) major,minor,particleclusternumber
major[d]=m_moments[2]*109/8
minor[d]=m_moments[3]*109/8
particleclusternumber[d]=l
d+=1
endfor
numberofparticles[l]=dimsize(m_moments,0)



//Display;AppendImage clu_NNimage
//ModifyImage clu_NNimage minRGB=NaN,maxRGB=0
//AppendToGraph W_BoundaryY vs W_BoundaryX
//ModifyGraph lsize=2,rgb=(0,0,0)
//ModifyGraph lsize=4
///ModifyGraph noLabel=2
//ModifyGraph tick=3,userticks=0
//ColorScale/C/N=text0/A=RC/E image=clu_NNimage
//ModifyImage clu_NNimage ctab= {1,*,Rainbow,1}
//ColorScale/C/N=text0/A=RC/E image=clu_NNimage
//string tosave1=path+"clust"+num2str(l)+".png"
//SavePICT/E=-5/B=576/o as tosave1



endfor

string tosave2=path+"All_Particle_Information_all.txt"
Save/O/J/M="\n"/W Major,Minor,Particleclusternumber  as tosave2

string tosave3=path+"Particles per Cluster_all.txt"
Save/O/J/M="\n"/W numberofparticles  as tosave3





  end
  
  
function checkparticle(num)
variable num
variable d
variable l

variable clustnumber=num
kill()
string xlist="x"+num2str(clustnumber)
string ylist="y"+num2str(clustnumber)
string dllist="dl"+num2str(clustnumber)
duplicate/o $xlist,tempx
duplicate/o $ylist,tempy
duplicate/o $dllist,tempd
variable scale=8
wavestats/q tempx
variable xmin=round(v_min)
variable xmax=round(v_max)
variable xdiff=round(v_max-v_min)

wavestats/q tempy
variable ymin=round(v_min)
variable ymax=round(v_max)
variable ydiff=round(v_max-v_min)

make/o/n=(scale*(xdiff+2),scale*(ydiff+2)) clu_DLimage=0,clu_Limage=0,clu_NNimage=0
variable a
for(a=0;a<(dimsize(tempx,0));a+=1)							// Populate the matrix
	variable xpos=round(scale*(tempx[a]-xmin+1))							// Get the co-ordinates
	variable ypos=round(scale*(tempy[a]-ymin+1))		
					
	clu_DLimage[xpos][ypos]=Clu_DLimage[xpos][ypos]+tempd[a]							// Populate the matrix
	clu_Limage[xpos][ypos]+=1
endfor

variable b

for(a=0;a<(dimsize(clu_dlimage,0));a+=1)
	for(b=0;b<(dimsize(clu_dlimage,1));b+=1)
		if(Clu_Limage[a][b]>0)
			Clu_NNimage[a][b]=Clu_DLimage[a][b]/clu_Limage[a][b]
		else
			clu_NNimage[a][b]=0
		endif
	endfor
endfor

ImageThreshold/T=30/I clu_nnimage

wave M_imagethresh

ImageAnalyzeParticles/E/W/Q/M=1/A=5 stats,M_imageThresh
newimage clu_nnimage



wave w_boundaryx,w_boundaryy
Display;AppendImage clu_NNimage
ModifyImage clu_NNimage minRGB=NaN,maxRGB=0
AppendToGraph W_BoundaryY vs W_BoundaryX
ModifyGraph lsize=2,rgb=(0,0,0)
ModifyGraph lsize=4
ModifyGraph noLabel=2
ModifyGraph tick=3,userticks=0
ColorScale/C/N=text0/A=RC/E image=clu_NNimage
ModifyImage clu_NNimage ctab= {1,*,Rainbow,1}
ColorScale/C/N=text0/A=RC/E image=clu_NNimage












  end
  
  
  
  /////////////// POST ANALYSIS FOR APTAMER PAPER////////////////
  
function post_analyse(first,last)
variable first				// First folder to look in
variable last				// Last folder to look in
setdatafolder root:
variable length=(last-first)
make/o/n=(length+1) number_of_clusters,pointsperclusterave,pointsperclustersdev,fold,resolution_ave,resolution_sd,aveNN,sdevNN,AveDist,SdevDist,Prec_ave,Prec_SD
make/o/n=(length+1) pointspercluster_med,resolution_med,NN_med,Dist_med,Prec_med
variable a,b,c

for(a=first;a<(last+1);a+=1)
	string tofolder=num2str(a)
	setdatafolder $tofolder
	 wave cluster 				// cluster address
	 wavestats/q cluster
	 variable numberofclust=v_max
	
	wave points_per_cluster
	wavestats/q points_per_cluster
	variable pointsave=v_avg
	variable pointssdev=v_sdev
	
	wave reff
	wavestats/q reff
	variable redav=v_avg
	variable refsd=v_sdev
	
	wave NumNN
	wavestats/q NumNN
	variable NNav=v_avg
	variable NNsd=v_sdev
	
	wave dl_distance
	duplicate/o DL_distance,NN_distance_ave
	
	
	
	for(c=0;c<(dimsize(DL_distance,0));c+=1)
		NN_distance_ave[c]=dl_distance[c]/5
	endfor
	
	
	wavestats/q NN_distance_ave
	variable DLav=v_avg
	variable DLsd=v_sdev
	
	
	wave Precision__nm_
	wavestats/q Precision__nm_
	variable precav=v_avg
	variable precsd=v_sdev
	
	
	variable medpoint=median(points_per_cluster)
	variable medres=median(reff)
	variable medNN=median(NN_distance_ave)
	variable meddist=median(dl_distance)
	variable medprec=median(precision__nm_)
	
	setdatafolder root:
	fold[b]=a
	number_of_clusters[b]=numberofclust
	pointsperclusterave[b]=pointsave
	pointsperclustersdev[b]=pointssdev
	resolution_ave[b]=redav
	resolution_sd[b]=refsd
	aveNN[b]=NNav
	sdevNN[b]=NNsd
	AveDist[b]=DLav
	SdevDist[b]=DLsd
	Prec_ave[b]=precav
	Prec_SD[b]=precsd
	pointspercluster_med[b]=medpoint
	resolution_med[b]=medres
	NN_med[b]=mednN
	Dist_med[b]=meddist
	Prec_med[b]=medprec




b+=1






endfor





















end





function post_analyse_clust(first,last)
variable first				// First folder to look in
variable last				// Last folder to look in
setdatafolder root:
variable length=(last-first)
make/o/n=(length+1) pointspercluster_c,pointspercluster_nc,resolution_c,resolution_nc,NN_c,NN_nc,AveDist_c,Avedist_nc,Prec_c,Prec_nc

variable a,b,c,d,e,f

for(a=first;a<(last+1);a+=1)
	string tofolder=num2str(a)
	setdatafolder $tofolder
	 wave cluster 				// cluster address
	 wavestats/q cluster
	 variable numberofclust=v_max
	wave coincident_clusters
	wave points_per_cluster
	wave numnn
	wave reff
		wave dl_distance
	duplicate/o DL_distance,NN_distance_ave
	
	
	
	for(c=0;c<(dimsize(DL_distance,0));c+=1)
		NN_distance_ave[c]=dl_distance[c]/5
	endfor
	
	
	make/o/n=1 points__c,points__nc,resolution__c,resolution__nc,NN__c,NN__nc,Dist__c,Dist__nc
		d=0
		e=0
		f=0
		for(c=0;c<(dimsize(coincident_clusters,0));c+=1)
			if(coincident_clusters[c]>0)
				redimension/n=(d+1) points__c,resolution__c,NN__c,Dist__c
				points__c[d]=points_per_cluster[c]
				resolution__c[d]=reff[c]
				NN__c[d]=numnn[c]
				dist__c[d]=NN_distance_ave[c]
				d+=1
			else
				redimension/n=(d+1) points__nc,resolution__nc,NN__nc,Dist__nc
				points__nc[e]=points_per_cluster[c]
				resolution__nc[e]=reff[c]
				NN__nc[e]=numnn[c]
				dist__nc[e]=NN_distance_ave[c]
				e+=1
			endif
		endfor
				
	
	
	
	wavestats/q points__c
	variable pointsave_c=v_avg
	
	wavestats/q resolution__c
	variable ref_c=v_sdev
	
	wavestats/q NN__c
	variable NN_c_=v_avg
	
	wavestats/q dist__c
	variable dist_c=v_avg
	
	wavestats/q points__nc
	variable pointsave_nc=v_avg
	
	wavestats/q resolution__nc
	variable ref_nc=v_sdev
	
	wavestats/q NN__nc
	variable NN_nc_=v_avg
	
	wavestats/q dist__nc
	variable dist_nc=v_avg
	
	
	setdatafolder root:
pointspercluster_c[b]=pointsave_c
pointspercluster_nc[b]=pointsave_nc
resolution_c[b]=ref_c
resolution_nc[b]=ref_nc
NN_c[b]=nn_c_
NN_nc[b]=nn_nc_
AveDist_c[b]=dist_c
Avedist_nc[b]=dist_nc





b+=1






endfor

end




function savedlplot(path)

string path
string tosave=path+"forrender.txt"
wave all_Frame, all_origX, all_origY, all_origValue, all_Error, all_Noise,all_Background,all_Signal, all_Angle, all_XW, all_YW, all_X_SD, all_Y_SD, all_Precision__nm_,numnn

make/o/n=(dimsize(all_frame,0)) used=1

//Save/J/M="\n"/O/W all_Frame, all_origX, all_origY, all_origValue, all_Error, all_Noise,all_Background,all_Signal, all_Angle, all_XW, all_YW, all_X_SD, all_Y_SD, all_Precision__nm_ as tosave
variable f1 
Open f1 as tosave


wfprintf f1, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t\r" all_Frame,all_xw,all_yw,all_X_SD,all_Y_SD,all_background,all_signal,numnn,used
Close f1





end









function DLclusterplot(path)			// This needs a prior code that outputs all localisations as xw, yw, precision. 
string path
wave all_xw,all_yw,all_precision__nm_,numNN
								// I only need the x and y-coords and precision values.
variable scale=8										// Change this for different scale.
variable size=scale*512								// This is the size of the matrix to make.
		// Make matrices for data. 

variable a,b,c											// Some originally named variables.
	

///// THIS PART FOR 2D GAUSS- width = precision image ///////////////////

variable initial=0

variable maximumnn=100
variable increment=round((maximumnn-initial)/100)
variable nnv
print increment


make/o/n=(11,11) tempgauss						// This is where I'm going to store each Gaussian that's fitted. Need to make an odd number in size to ensure there is a central pixel. 
make/o/n=(size,size) tempimagegauss=0
newimage tempimagegauss
ModifyGraph tick=3,noLabel=2
SetDrawLayer ProgFront
SetDrawEnv fillfgc= (0,65535,0),linethick=0,fillpat=1,xcoord= top,ycoord= left, save
ModifyImage tempimagegauss ctab= {*,1,Grays,0}
ColorTab2Wave rainbow
c=0
wave m_colors
duplicate/o m_colors,colorwave
for(nnv=initial;nnv<maximumnn;nnv+=increment)



	duplicate/o all_precision__nm_,precisionpixel			// Notekeeping- I did this just so that I could confirm that the precisions are correct when converted to pixels (from nm)
	for(a=0;a<(dimsize(all_xw,0));a+=1)
		if(numNN[a]<=(nnv+increment) && numNN[a]>nnv)
			variable xcoord=round(scale*(all_xw[a]))
			variable ycoord=round(scale*(all_yw[a]))
			variable precis=all_precision__nm_[a]
		// Need to convert precision to pixels
			variable pixelsize=109						// Change pixel size here if different.
			variable pixelsizeinexpanded=109/scale	// This is new pixel size with different scale. 
			variable precispixel=precis/pixelsizeinexpanded			// This is now the precision in pixels rather than nm. 

			//precisionpixel[a]=precispixel								// Store for notekeeping. 
	
			//duplicate/o tempgauss,tempfit								// I need to make a wave to store the fit in. Idiosynchratic nature of Igor!!

			//K0 = 0;K1 = 1;K2 = 5;K3 = precispixel;K4 = 5;K5 = precispixel;K6 = 0;				// These are the fixed variables to fit the 2D Gaussian to - amplitude = 1, x-width,ywidth = precision (in pixels)

			//CurveFit/H="1111111"/q Gauss2D tempgauss /D=tempfit 				// Fit the curve, and output the fit to tempfit wave. 

			//for(b=0;b<11;b+=1)							// Place this in the original gaussimage.
				//f//or(c=0;c<11;c+=1)
					//tempimagegauss[xcoord-5+b][ycoord-5+c]=tempimagegauss[xcoord-5+b][ycoord-5+c]+tempfit[b][c]				// Note that it's an add function - allows gauss to overlap. 
				//endfor
			//endfor
			variable r=colorwave[99-nnv][0]
			variable g=colorwave[99-nnv][1]
			variable bl=colorwave[99-nnv][2]
			
			SetDrawEnv fillfgc= (r,g,bl),linethick=0,fillpat=1,xcoord= top,ycoord= left, save
			drawoval	(xcoord-precispixel),(ycoord-precispixel),(xcoord+precispixel),(ycoord+precispixel)
			
endif
endfor
c+=1
endfor



end

function circs()
wave tempimagegauss
newimage tempimagegauss
SetDrawLayer ProgFront
SetDrawEnv fillfgc= (0,65535,0)
DrawRect 0,50,132,500


ModifyImage tempimagegauss ctab= {9000,20000,Grays,0}
ModifyGraph tick=3,noLabel=2
SetDrawLayer ProgFront
SetDrawEnv linefgc= (0,65535,0),fillpat= 1,xcoord= top,ycoord= left, save
SetDrawEnv fillfgc= (0,65535,0)
DrawRect 0,50,132,500
end




function plotallSR(path)
string path
variable clustnumber
wave NN_image
NewImage/K=0 NN_Image
ModifyImage NN_Image ctab={*,*,Rainbow,0}
ModifyImage NN_Image ctab= {*,*,Rainbow,1}
SetDrawLayer ProgFront
SetDrawEnv linefgc= (65535,65535,65535),fillpat= 0,xcoord= top,ycoord= left, save

variable c
wave xw
wavestats xw
variable maxi=v_max
for(c=0;c<(maxi);c+=1)
string xlist="x"+num2str(c)
string ylist="y"+num2str(c)
string dllist="dl"+num2str(c)
duplicate/o $xlist,tempx
duplicate/o $ylist,tempy
wavestats/q tempx
//variable xmin=round(v_min-100)*8
//variable xmax=round(v_max-100)*8
//variable xdiff=round(v_max-v_min)

variable xmin=round(v_min)*8
variable xmax=round(v_max)*8
variable xdiff=round(v_max-v_min)

wavestats/q tempy
//variable ymin=round(v_min-100)*8
//variable ymax=round(v_max-100)*8
//variable ydiff=round(v_max-v_min)

variable ymin=round(v_min)*8
variable ymax=round(v_max)*8
variable ydiff=round(v_max-v_min)
wave dl_distance,reff
DrawRect xmin,ymin,xmax,ymax


string dist=num2str(dl_distance[c]/5)
string res=num2str(reff[c])

string clu=num2str(c)
SetDrawEnv textrgb= (65535,65535,65535),fsize= 10;DelayUpdate
DrawText xmax,ymax,clu
//print xmin
//print xmax
endfor


end





function convert(yco,xco)

variable xco,yco
kill()
  combineall()
  clusterplotpaperall()
// First step - flip horizontal

variable xco1=(512*8)-xco
variable yco1=yco


//variable yco2=(8*512-xco1)
//variable xco2=(yco1)
variable yco2=(yco)
 variable xco2=(xco)
print yco2
print xco2

SetAxis/R bottom ((xco2-40)/8),((xco2+40)/8);DelayUpdate
SetAxis left ((yco2-40)/8),((yco2+40)/8)
modifyGraph msize=4,useMrkStrokeRGB=1
ModifyGraph msize=4,useMrkStrokeRGB=1

ColorScale/K/N=text1
ModifyGraph tick=3,minor(bottom)=0,noLabel=2
ModifyGraph width=100,height=100
ModifyGraph msize=2
ModifyGraph useMrkStrokeRGB=0
ModifyGraph msize(yco27)=1
ModifyGraph mrkThick('yco34.5')=0.1,useMrkStrokeRGB('yco34.5')=1
ModifyGraph msize('yco46.5')=1
ModifyGraph msize=1,mrkThick=0.1
ModifyGraph msize=1.5

ModifyGraph useMrkStrokeRGB=1
SavePICT/E=-2
plotall()

SetAxis top ((xco2-50)),((xco2+50));DelayUpdate
SetAxis left ((yco2-50)),((yco2+50))


end



function drawboxes()
wave imagegauss
make/o/n=6 xco={1588,1657,1736,2265,2257,2353},yco={2061,2415,2487,2952,1882,1841}

newimage imagegauss
ModifyImage imagegauss ctab= {0,50,Rainbow,0}
ModifyGraph tick=3,noLabel=2
ModifyImage imagegauss ctab= {0,50,YellowHot,0}



SetDrawLayer ProgFront
SetDrawEnv linethick=1,fillpat=0,xcoord= top,ycoord=left,linefgc=(65535,65535,65535), save


variable a

for(a=0;a<(dimsize(xco,0));a+=1)
variable left=xco[a]-40
variable top=yco[a]-40
variable right=xco[a]+40
variable bottom=yco[a]+40
SetDrawLayer ProgFront
DrawRect left,top,right,bottom
endfor
end

