#pragma rtGlobals=3		// Use modern global access method and strict wave access.

function multiple()
variable numberoffolders=3
make/o/n=(numberoffolders)/T filelist
filelist[0]="RAID:13_02_Cells:20170213_AptPaint_Cell_GS:561_1:0:"
filelist[1]="RAID:13_02_Cells:20170213_AptPaint_Cell_GS:561_2:0:"
filelist[2]="RAID:13_02_Cells:20170213_AptPaint_Cell_GS:561_4:0:"





variable f

for(f=0;f<(numberoffolders);f+=1)
	setdatafolder root:
	string folder=num2str(f)
	PRINT FOLDER
	newdatafolder/s $folder
	
	string path=filelist[f]
	
	loadmult(path)
	imageclusterplot(path)	
	
	clusteranal()
	setdatafolder root:
endfor



end

//////////LOAD MULTIPLE FILES////////////


function loadmult(path)
string path
string load1=path+"cluster.txt"
print load1
LoadWave/J/D/W/K=0/A load1

wave xw
duplicate/o xw,cluster
killwaves xw

string load2=path+"fitresults.txt"
print load2

LoadWave/J/D/W/K=0/A load2
end


///////////Output images//////////////

function imageclusterplot(path)			// This needs a prior code that outputs all localisations as xw, yw, precision. 
string path
wave xw,yw,precision__nm_,cluster
duplicate/o precision__nm_,precision									// I only need the x and y-coords and precision values.
variable scale=8										// Change this for different scale.
variable size=scale*512								// This is the size of the matrix to make.
make/o/n=(size,size) Image=0,imagegauss=0		// Make matrices for data. 

variable a,b,c											// Some originally named variables.
	
/////////// THIS PART TO PLOT LOCALISATIONS ONLY /////////////////////////	
		
for(a=0;a<(dimsize(xw,0));a+=1)					// Go through all of localisations.
	if(cluster[a]>0)
		variable xcoord=round(scale*(xw[a]))		// Multiply by 8 to fit in matrix. 
		variable ycoord=round(scale*(yw[a]))		
		image[xcoord][ycoord]+=1						// Add 1 to the matrix- localisation only. 
	endif
endfor


///// THIS PART FOR 2D GAUSS- width = precision image ///////////////////

make/o/n=(11,11) tempgauss						// This is where I'm going to store each Gaussian that's fitted. Need to make an odd number in size to ensure there is a central pixel. 

duplicate/o precision,precisionpixel			// Notekeeping- I did this just so that I could confirm that the precisions are correct when converted to pixels (from nm)
for(a=0;a<(dimsize(xw,0));a+=1)
if(cluster[a]>0)
	xcoord=round(scale*(xw[a]))
	ycoord=round(scale*(yw[a]))
	variable precis=precision[a]
		// Need to convert precision to pixels
	variable pixelsize=109						// Change pixel size here if different.
	variable pixelsizeinexpanded=109/scale	// This is new pixel size with different scale. 
	variable precispixel=precis/pixelsizeinexpanded			// This is now the precision in pixels rather than nm. 

	precisionpixel[a]=precispixel								// Store for notekeeping. 
	
	duplicate/o tempgauss,tempfit								// I need to make a wave to store the fit in. Idiosynchratic nature of Igor!!

	K0 = 0;K1 = 1;K2 = 5;K3 = precispixel;K4 = 5;K5 = precispixel;K6 = 0;				// These are the fixed variables to fit the 2D Gaussian to - amplitude = 1, x-width,ywidth = precision (in pixels)

	CurveFit/H="1111111"/q Gauss2D tempgauss /D=tempfit 				// Fit the curve, and output the fit to tempfit wave. 

	for(b=0;b<11;b+=1)							// Place this in the original gaussimage.
		for(c=0;c<11;c+=1)
			imagegauss[xcoord-5+b][ycoord-5+c]=imagegauss[xcoord-5+b][ycoord-5+c]+tempfit[b][c]				// Note that it's an add function - allows gauss to overlap. 
		endfor
	endfor
endif
endfor

string tosave=path+"imagegauss2.txt"

Save/J/M="\n" imagegauss as tosave // Saves the file. 

end


///////////////////////////////////////




function load_dan()
LoadWave/J/D/W/K=0/A

wave wave0,wave1,wave2

duplicate/o wave0,xw
duplicate/o wave1,yw
duplicate/o wave2,precision__nm_

end




////// Without clusters

function imageplot()			// This needs a prior code that outputs all localisations as xw, yw, precision. 
wave xw,yw,precision__nm_
duplicate/o precision__nm_,precision									// I only need the x and y-coords and precision values.
variable scale=8										// Change this for different scale.
variable size=scale*512								// This is the size of the matrix to make.
make/o/n=(size,size) Image=0,imagegauss=0		// Make matrices for data. 

variable a,b,c											// Some originally named variables.
	
/////////// THIS PART TO PLOT LOCALISATIONS ONLY /////////////////////////	
		
for(a=0;a<(dimsize(xw,0));a+=1)					// Go through all of localisations.
		variable xcoord=round(scale*(xw[a]))		// Multiply by 8 to fit in matrix. 
		variable ycoord=round(scale*(yw[a]))		
		image[xcoord][ycoord]+=1						// Add 1 to the matrix- localisation only. 
endfor


///// THIS PART FOR 2D GAUSS- width = precision image ///////////////////

make/o/n=(11,11) tempgauss						// This is where I'm going to store each Gaussian that's fitted. Need to make an odd number in size to ensure there is a central pixel. 

duplicate/o precision,precisionpixel			// Notekeeping- I did this just so that I could confirm that the precisions are correct when converted to pixels (from nm)
for(a=0;a<(dimsize(xw,0));a+=1)
	xcoord=round(scale*(xw[a]))
	ycoord=round(scale*(yw[a]))
	variable precis=precision[a]
		// Need to convert precision to pixels
	variable pixelsize=109						// Change pixel size here if different.
	variable pixelsizeinexpanded=109/scale	// This is new pixel size with different scale. 
	variable precispixel=precis/pixelsizeinexpanded			// This is now the precision in pixels rather than nm. 

	precisionpixel[a]=precispixel								// Store for notekeeping. 
	
	duplicate/o tempgauss,tempfit								// I need to make a wave to store the fit in. Idiosynchratic nature of Igor!!

	K0 = 0;K1 = 1;K2 = 5;K3 = precispixel;K4 = 5;K5 = precispixel;K6 = 0;				// These are the fixed variables to fit the 2D Gaussian to - amplitude = 1, x-width,ywidth = precision (in pixels)

	CurveFit/H="1111111"/q Gauss2D tempgauss /D=tempfit 				// Fit the curve, and output the fit to tempfit wave. 

	for(b=0;b<11;b+=1)							// Place this in the original gaussimage.
		for(c=0;c<11;c+=1)
			imagegauss[xcoord-5+b][ycoord-5+c]=imagegauss[xcoord-5+b][ycoord-5+c]+tempfit[b][c]				// Note that it's an add function - allows gauss to overlap. 
		endfor
	endfor

endfor

Save/J/M="\n" imagegauss // Saves the file. 

end


function clusteranal()
wave cluster,xw,yw,precision

variable cl,rows

wavestats cluster
variable length=v_max
make/o/n=(length) numberofpoints
make/o/n=1 temp_prec,temp_x,temp_y
for(cl=0;cl<length;cl+=1)
variable count=0
	for(rows=0;rows<(dimsize(cluster,0));rows+=1)
		if(cluster[rows]==(cl+1))
			
			
			redimension/n=(count+1) temp_prec,temp_x,temp_y
			temp_prec[count]=precision[rows]
			temp_x[count]=xw[rows]
			temp_y[count]=yw[rows]
					
			
			count+=1
		endif
	endfor
string name1=num2str(cl)+"xw"
string name2=num2str(cl)+"yw"
string name3=num2str(cl)+"precision"
numberofpoints[cl]=count
duplicate/o temp_prec,$name3
duplicate/o temp_x,$name1
duplicate/o temp_y,$name2





endfor

Make/N=25/O numberofpoints_Hist;DelayUpdate
Histogram/B={0,25,25} numberofpoints,numberofpoints_Hist;DelayUpdate
Display numberofpoints_Hist
SetAxis bottom *,1000
ModifyGraph mode=5
print median(numberofpoints)
end


function clusterloaddan()
LoadWave/J/D/W/K=0/A

wave xw
duplicate/o xw,cluster
killwaves xw

LoadWave/J/D/W/K=0/A
wave wave0,wave1,wave2

duplicate/o wave0,xw
duplicate/o wave1,yw
duplicate/o wave2,precision__nm_


end





function plot_number(num)
variable num

string name1=num2str(num)+"xw"
string name2=num2str(num)+"yw"
string name3=num2str(num)+"precision"

duplicate/o $name3,temp_prec
duplicate/o $name1,temp_x
duplicate/o $name2,temp_y


make/o/n=(11,11) tempgauss						// This is where I'm going to store each Gaussian that's fitted. Need to make an odd number in size to ensure there is a central pixel. 

variable a,b,c,scale=8
make/o/n=(scale*512,scale*512) t_imagegauss=0
duplicate/o temp_prec,precisionpixel			// Notekeeping- I did this just so that I could confirm that the precisions are correct when converted to pixels (from nm)
for(a=0;a<(dimsize(temp_x,0));a+=1)

	variable xcoord=round(scale*(temp_x[a]))
	variable ycoord=round(scale*(temp_y[a]))
	variable precis=temp_prec[a]
		// Need to convert precision to pixels
	variable pixelsize=109						// Change pixel size here if different.
	variable pixelsizeinexpanded=109/scale	// This is new pixel size with different scale. 
	variable precispixel=precis/pixelsizeinexpanded			// This is now the precision in pixels rather than nm. 

	precisionpixel[a]=precispixel								// Store for notekeeping. 
	
	duplicate/o tempgauss,tempfit								// I need to make a wave to store the fit in. Idiosynchratic nature of Igor!!

	K0 = 0;K1 = 1;K2 = 5;K3 = precispixel;K4 = 5;K5 = precispixel;K6 = 0;				// These are the fixed variables to fit the 2D Gaussian to - amplitude = 1, x-width,ywidth = precision (in pixels)

	CurveFit/H="1111111"/q Gauss2D tempgauss /D=tempfit 				// Fit the curve, and output the fit to tempfit wave. 

	for(b=0;b<11;b+=1)							// Place this in the original gaussimage.
		for(c=0;c<11;c+=1)
			t_imagegauss[xcoord-5+b][ycoord-5+c]=t_imagegauss[xcoord-5+b][ycoord-5+c]+tempfit[b][c]				// Note that it's an add function - allows gauss to overlap. 
		endfor
	endfor

endfor

wavestats temp_x

variable xmax=8*v_max
variable xmin=8*v_min
variable xrange=(xmax-xmin)
wavestats temp_y

variable ymax=8*v_max
variable ymin=8*v_min
variable yrange=(xmax-xmin)

if(yrange>xrange)
variable length=yrange
else
length=xrange
endif



newimage t_imagegauss
//AppendImage t_imagegauss
ModifyGraph tick=3,noLabel=2
SetAxis/R top (xmin-40),(xmin+length+40)
SetAxis/R left (ymin-40),(ymin+length+40)
ModifyImage t_imagegauss ctab= {*,*,YellowHot,0}
ModifyGraph width=283.465,height=283.465
end