run("Set Measurements...", "area center fit feret's redirect=None decimal=5");
for(img_num=4118; img_num<=4119; img_num++){
	folder = ".";
	if(img_num<10)
		img_name = "000"+img_num;
	else if(img_num<100)
		img_name = "00"+img_num;
	else if(img_num<1000)
		img_name = "0"+img_num;
	else
		img_name = img_num;
    file = folder+"/IMG_"+img_name+".JPG";
	open(file);
	run("Unsharp Mask...", "radius=4 mask=0.90");
	min=newArray(3);
	max=newArray(3);
	filter=newArray(3);
	a=getTitle();
	run("HSB Stack");
	run("Convert Stack to Images");
	selectWindow("Hue");
	rename("0");
	selectWindow("Saturation");
	rename("1");
	selectWindow("Brightness");
	rename("2");
	min[0]=0;
	max[0]=255;
	filter[0]="pass";
	min[1]=0;
	max[1]=255;
	filter[1]="pass";
	min[2]=0;
	max[2]=250;
	filter[2]="pass";
	for (i=0;i<3;i++){
	  selectWindow(""+i);
	  setThreshold(min[i], max[i]);
	  run("Convert to Mask");
	  if (filter[i]=="stop")  run("Invert");
	}
	imageCalculator("AND create", "0","1");
	imageCalculator("AND create", "Result of 0","2");
	for (i=0;i<3;i++){
	  selectWindow(""+i);
	  close();
	}
	selectWindow("Result of 0");
	close();
	selectWindow("Result of Result of 0");
	rename(a);
    //makeOval(530, 140, 1596, 1596);
    makeOval(270, 100, 1650, 1650);
	run("Make Inverse");
	run("Fill", "slice");
	run("Invert");
	run("Make Binary");
	selectWindow("IMG_"+img_name+".JPG");
	saveAs("Jpeg", folder+"/binary_"+img_name+".jpg");
	run("Analyze Particles...", "size=100-1500 circularity=0-10 show=Outlines display clear summarize record");
	selectWindow("Drawing of binary_"+img_name+".jpg");
	saveAs("Jpeg", folder+"/Ellipses_"+img_name+".jpg");
	run("Input/Output...", "jpeg=85 gif=-1 file=.txt use_file copy_row save_row");
	saveAs("Results", folder+"/rods_"+img_name);
	run("Close All");
}
