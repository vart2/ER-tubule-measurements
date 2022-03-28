
scale = 3.5;			//scale in nm/pixel
stackName = "mut3 w150"	//name of image stack
resultsFolder = "/Users/vishal/Desktop/Project/mut3 w150/2.5 blur, global, minferet"	//Path to folder to save results
exclusionRoi = "/Users/vishal/Desktop/Project/mut3 w150/0001-1566-3735.roi"	//Path to ROI for exclusion of glia


feretMaxLim = 125;	//max diameter value in nm
feretMinLim = 10;	//min diameter value in nm
disLim = 50;		//distance tubules can move between slices in nm

feretMaxLim = feretMaxLim/scale;
feretMinLim = feretMinLim/scale;
disLim = disLim/scale;

f = 0;
selectionRadius = disLim + feretMaxLim;
selectionDiameter = 2 * selectionRadius;
disLim = Math.pow(disLim, 2);
setOption("ExpandableArrays", true);
var rois = newArray;
var distance = newArray;
var xx = newArray;
var yy = newArray;

function checkRois (array, value) {
	check = false;
	if (array.length!=0) {
		for (i = 0; i < array.length; i++) {
			if (array[i]==value) {
				roiManager("deselect");
				roiManager("delete");
				check = true;
				break;
			}
		}
	}
	if (check == true) {return true;}
	rois[rois.length] = name;
	return false;
}

function tubuleDistance() { 
	run("Select None");
	xx = Array.trim(xx, 0);
	yy = Array.trim(yy, 0);
	distance = Array.trim(distance, 0);
	makeOval(xCoord-selectionRadius, yCoord-selectionRadius, selectionDiameter, selectionDiameter);
	run("Find Maxima...", "prominence=100 output=[Point Selection]");
	getSelectionCoordinates(xx, yy);
	if (xx.length==0) {return disLim+1;}
	for (i = 0; i < xx.length; i++) {
		distance[i] = Math.pow(xx[i]-xCoord, 2) + Math.pow(yy[i]-yCoord, 2);
	}
	Array.getStatistics(distance, min, max, mean, stdDev);
	return min;
}

function selectTubule() {
	minima = Array.findMinima(distance, 0);
	if (minima.length==0) {doWand(xx[0], yy[0], 0.0, "4-connected");}
	else {doWand(xx[minima[0]], yy[minima[0]], 0.0, "4-connected");}
	roiManager("add");
	roiManager("select", roiManager("count")-1);
}

function measure(startSlice, tubuleNumber) {
	if (roiManager("count")==0) {}
	else {
		if (roiManager("count")>2) {
		roiManager("deselect");
		roiManager("Save", resultsFolder+"/"+stackName+" RoiSet"+startSlice+", "+tubuleNumber+".zip");
		roiManager("measure");
		saveAs("Results", resultsFolder+"/"+stackName+" Results"+startSlice+", "+tubuleNumber+".csv");
		run("Clear Results");
		roiManager("delete");
		}
		else {
		roiManager("deselect");
		roiManager("delete");
		}
	}
}


for (k = 1; k < nSlices-1; k++) {
run("Select None");
open(exclusionRoi);
setSlice(k);
Roi.setPosition(k);
roiManager("add");
roiManager("Select", 0);
run("Find Maxima...", "prominence=100 output=[Point Selection]");
getSelectionCoordinates(xpoints, ypoints);
roiManager("delete");
for (j = 0; j < xpoints.length; j++) {
setSlice(k);
doWand(xpoints[j], ypoints[j], 0.0, "4-connected");
roiManager("add");
roiManager("select", roiManager("count")-1);
f = getValue("MinFeret");

while ((f<feretMaxLim)&&(f>feretMinLim)) {
roiName = Roi.getName();
xCoord = substring(roiName, 10);
yCoord = substring(roiName, 5, 9);
if (getSliceNumber()==nSlices) {break;}
setSlice(getSliceNumber()+1);
if (tubuleDistance()>disLim) {
	if (getSliceNumber()==nSlices) {break;}
	setSlice(getSliceNumber()+1);
	if (tubuleDistance()>disLim) {break;}
	selectTubule();
	name = Roi.getName;
	if (checkRois(rois, name)) {break;}
	f = getValue("MinFeret");
	if ((f>=feretMaxLim)||(f<=feretMinLim)) {roiManager("delete"); break;}
}
else {
	selectTubule();
	name = Roi.getName;
	if (checkRois(rois, name)) {break;}
	f = getValue("MinFeret");
	if ((f<feretMaxLim)&&(f>feretMinLim)) {}
	else {
		roiManager("delete");
		if (getSliceNumber()==nSlices) {break;}
		setSlice(getSliceNumber()+1);
		if (tubuleDistance()>disLim) {break;}
		selectTubule();
		name = Roi.getName;
		if (checkRois(rois, name)) {break;}
		f = getValue("MinFeret");
		if ((f>=feretMaxLim)||(f<=feretMinLim)) {roiManager("delete"); break;}
	}
}
}

measure(k, j);
}
}
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
print(hour+":"+minute+":"+second);