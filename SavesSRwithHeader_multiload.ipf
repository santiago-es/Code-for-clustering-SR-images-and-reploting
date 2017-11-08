#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3		// Use modern global access method and strict wave access.


macro all()

multiple()

end

function multiple()
variable numberoffolders=14
variable numberoffiles=10
make/o/n=(numberoffolders)/T filelist

           filelist[0]="Birdbox III:Birdbox1_20160217:20161109_MitoAsyn:allprimaries:1:561:"

            filelist[1]="Birdbox III:Birdbox1_20160217:20161109_MitoAsyn:allprimaries:1:641:"

            filelist[2]="Birdbox III:Birdbox1_20160217:20161109_MitoAsyn:allprimaries:2:561:"

            filelist[3]="Birdbox III:Birdbox1_20160217:20161109_MitoAsyn:allprimaries:2:641:"

            filelist[4]="Birdbox III:Birdbox1_20160217:20161109_MitoAsyn:allprimaries:3:561:"

            filelist[5]="Birdbox III:Birdbox1_20160217:20161109_MitoAsyn:allprimaries:3:641:"

            filelist[6]="Birdbox III:Birdbox1_20160217:20161109_MitoAsyn:Monomerplusprimaries:2016-11-09_16-37-38:561:"

            filelist[7]="Birdbox III:Birdbox1_20160217:20161109_MitoAsyn:Monomerplusprimaries:2016-11-09_16-37-38:641:"

            filelist[8]="Birdbox III:Birdbox1_20160217:20161109_MitoAsyn:noprimaries:1:561:"

            filelist[9]="Birdbox III:Birdbox1_20160217:20161109_MitoAsyn:noprimaries:1:641:"

            filelist[10]="Birdbox III:Birdbox1_20160217:20161109_MitoAsyn:noprimaries:2:561:"

            filelist[11]="Birdbox III:Birdbox1_20160217:20161109_MitoAsyn:noprimaries:2:641:"

            filelist[12]="Birdbox III:Birdbox1_20160217:20161109_MitoAsyn:noprimaries:3:561:"

            filelist[13]="Birdbox III:Birdbox1_20160217:20161109_MitoAsyn:noprimaries:3:641:"







variable f

for(f=0;f<(numberoffolders);f+=1)
	setdatafolder root:
	string folder=num2str(f)
	PRINT FOLDER
	newdatafolder/s $folder
	
	string path=filelist[f]
	
	load(path,numberoffiles)
	
	setdatafolder root:
endfor



end





function load(path,numberoffolders)
string path
variable numberoffolders
string folder=path
print folder

string filename="fitresults.txt"

variable a,b,c
variable correct=0
//Load the files
make/o/n=1 all_Frame, all_origX, all_origY, all_origValue, all_Error, all_Noise, all_SNR,all_Background,all_Signal, all_Angle, all_XW, all_YW, all_X_SD, all_Y_SD, all_Precision__nm_
for(a=0;a<(numberoffolders);a+=1)

// If zeros in name:
if(a<10)
	string filetoload=folder+"0"+num2str(a)+":"+filename
	//string filetoload=folder+num2str(a)+":"+filename
else
	filetoload=folder+num2str(a)+":"+filename
endif
	LoadWave/q/o/J/D/W/A/K=0 filetoload
	print filetoload
	wave Source, Frame, origX, origY, origValue, Error, Noise, SNR,Background, Signal, Angle, XW, YW, X_SD,Y_SD, Precision__nm_
	
	for(b=0;b<(dimsize(frame,0));b+=1)
		if(precision__nm_[b]<20)
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
	endif

	endfor
		killwaves source
	killwaves frame
	killwaves origx
	killwaves origy
	killwaves origvalue
	killwaves error
	killwaves noise
	killwaves snr
	killwaves background
	killwaves signal
	killwaves angle
	killwaves xw
	killwaves yw
	killwaves x_sd
	killwaves y_sd
	killwaves precision__nm_
	wavestats/q all_frame
	correct=v_max


endfor

string tosave=folder+"allfitswithheader.txt"



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
fprintf f1, "#Frame\torigX\torigY\torigValue\tError\tNoise\tBackground\tSignal\tAngle\tX\tY\tXSD\tYSD\tPrecision\r"
wfprintf f1, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\r" all_Frame, all_origX, all_origY, all_origValue, all_Error, all_Noise,all_Background,all_Signal, all_Angle, all_XW, all_YW, all_X_SD, all_Y_SD, all_Precision__nm_
Close f1

string tosave2="Macintosh HD:Users:Mathew:Desktop:indexlist.txt"
variable f2
open/A f2 as tosave2
fprintf f2, "path[]=\"%s\"\r" folder
close f2


string tosave3=folder+"FitResults.txt"
variable f3
Open f3 as tosave3
fprintf f3, "#Frame\torigX\torigY\torigValue\tError\tNoise\tBackground\tSignal\tAngle\tX\tY\tXSD\tYSD\tPrecision\r"
wfprintf f3, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\r" all_Frame, all_origX, all_origY, all_origValue, all_Error, all_Noise,all_Background,all_Signal, all_Angle, all_XW, all_YW, all_X_SD, all_Y_SD, all_Precision__nm_
Close f3

end


