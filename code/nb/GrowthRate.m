(* ::Package:: *)

(* ::Section:: *)
(*Growth rate*)


(* ::Subsubsection:: *)
(*Import data*)


(* ::Input:: *)
(*tabqx=Import["../data/Dump_Growth_Rate_Precession.hf5",{"Datasets","tabqx"}];*)
(*tabGR=Import["../data/Dump_Growth_Rate_Precession.hf5",{"Datasets","tabGR"}];*)
(**)
(*GRTable=Partition[{tabqx,tabGR}//Transpose//Flatten,3];*)
(*a=Log10[GRTable[[All,1]]];*)
(*b=GRTable[[All,2]];*)
(*GRTable[[All,2]]=a;*)
(*GRTable[[All,1]]=b;*)


(* ::Subsubsection:: *)
(*Plot data*)


(* ::Input:: *)
(*(*Generic options for the plots.*)*)
(*fontname = "Helvetica";*)
(*fontsize = 15;*)
(*basestyleplot = {FontFamily -> fontname,FontSize->fontsize};*)
(*labelstyleplot = {FontFamily -> fontname, FontSize -> fontsize};*)
(*imagesize = Large;(*350;*)*)
(*(*Making the plot.*)*)
(*(*Range of the plot.*)*)
(*qmin = 0.01;*)
(*qmax =1;*)
(*xmin = 0.0;*)
(*xmax = 0.4;*)
(*zmin =0;*)
(*zmax = 2;*)
(*(*Axes of the plot.*)*)
(**)
(*ax = Row[{Style["log(q)",fontname,fontsize]}];*)
(*ay =  Row[{Style["x",fontname,fontsize]}];*)
(*(*.*)*)
(*axstyle = {Black,Thickness[0.003]};*)
(*tickstyle = {Black,Thickness[0.003]};*)
(*(*.*)*)
(*qticks =  {{-3,"-3.0"},{-2.75,""},{-2.5,"-2.5"},{-2.25,""},{-2,"-2.0"},{-1.75,""},{-1.5,"-1.5"},{-1.25,""},{-1,"-1.0"},{-0.75,""},{-0.5,"-0.5"},{-0.25,""},{0,"0.0"}};*)
(*xticks =  {{0.0,"0.0"},{0.1,"0.1"},{0.2,"0.2"},{0.3,"0.3"},{0.4,"0.4"},{0.5,"0.5"},{0.6,"0.6"},{0.7,"0.7"},{0.8,"0.8"},{0.9,"0.9"},{1.0,"1.0"}};*)
(*GRTicks={{0,"0"},{0.5,""},{0.25,""},{0.5,"0.5"},{0.75,""},{1,"1"},{1.25,""},{1.5,"1.5"},{1.75,""},{2,"2"}};*)


(* ::Input:: *)
(*col0=Join[ColorData["SunsetColors"]/@(Range[0,15]/15),{White}]//Reverse;*)
(*nn=Length[col0]-1;*)
(*col[x_]=Blend[Transpose[{Range[0,nn]/nn,col0}],x];*)
(*img=Import["spiral.png"];*)
(*im2=ColorNegate[img];*)
(*img2= Import["spiral2.png"];*)


(* ::Input:: *)
(*contours=Table[0.1+i*0.1,{i,0,19}];*)
(*pGR=ListContourPlot[GRTable,Contours->contours,PlotRange->{{xmin,xmax},Log10[{qmin,qmax}],{zmin,zmax}},ContourStyle->None,ColorFunction->col,FrameLabel->{ay,ax},BaseStyle->basestyleplot,LabelStyle->labelstyleplot,PlotLegends->Placed[BarLegend[{col,{zmin,zmax}},LegendLabel->"Growth rate \!\(\*SqrtBox[TagBox[FractionBox[SuperscriptBox[SubscriptBox[\"a\", \"d\"], \"3\"], *)
(*RowBox[{\"G\", \" \", \"M\"}]],*)
(*DisplayForm]]\)\!\(\*SubscriptBox[\(\[Omega]\), \(I\)]\) for \!\(\*TagBox[FractionBox[SubscriptBox[\(a\), \(d\)], SubscriptBox[\(a\), \(b\)]],*)
(*DisplayForm]\)=20, \!\(\*TagBox[FractionBox[SubscriptBox[\(a\), \(h\)], SubscriptBox[\(a\), \(d\)]],*)
(*DisplayForm]\)=2.8 and z=1.3",*)
(*Method->{FrameStyle->axstyle,TicksStyle->tickstyle,Ticks->GRTicks}],Above],ImageSize->imagesize,FrameTicks->{xticks,qticks},FrameStyle->{axstyle,axstyle,axstyle,axstyle},PlotRangePadding->None];*)
(*liDatag={{0.005,-1.7},{0.2,-1.0}};(*z=1.3*)*)
(*pDatag=ListPlot[liDatag,PlotRange->{{xmin,xmax},Log10[{qmin,qmax}]},PlotMarkers->{img2,Scaled[0.08]},PlotStyle->ColorData[11][5]];*)
(*liDatag2={{0.24,-0.34}};(*z=1.3*)*)
(*pDatag2=ListPlot[liDatag2,PlotRange->{{xmin,xmax},Log10[{qmin,qmax}]},PlotMarkers->{img,Scaled[0.08]},PlotStyle->ColorData[11][5]];*)
(*pGRData=Show[pGR,pDatag,pDatag2];*)
(*pGRData=Rasterize[pGRData,RasterSize->1000]*)
(*pDensGR=ListPlot3D[GRTable,PlotRange->{{xmin,xmax},Log10[{qmin,qmax}],{zmin,zmax}},ColorFunction->col,BaseStyle->basestyleplot,LabelStyle->labelstyleplot,ImageSize->imagesize,PlotRangePadding->None,AxesLabel->{ay,ax,"\!\(\*SqrtBox[TagBox[FractionBox[SuperscriptBox[\"a\", \"3\"], *)
(*RowBox[{\"G\", \" \", \"M\"}]],*)
(*DisplayForm]]\)\!\(\*SubscriptBox[\(\[Omega]\), \(I\)]\)"},Ticks->{xticks,qticks,GRTicks},MeshFunctions->{#3&}];*)
(*liData3D={{0,-1.7,0.05},{0.2,-1.0,0.17}};*)
(*pData3D=ListPointPlot3D[liData3D,PlotStyle->{Yellow,AbsolutePointSize[15]}];*)
(*liData3Db={{0.24,-0.34,0.33}};*)
(*pData3Db=ListPointPlot3D[liData3Db,PlotStyle->{Red,AbsolutePointSize[15]}];*)
(*pDensGRData=Show[pDensGR,pData3D,pData3Db];*)
(*pDensGRData=Rasterize[pDensGRData,RasterSize->1000]*)


(* ::Subsubsection:: *)
(*Save data*)


(* ::Input:: *)
(*Export["../graphs/GrowthRate.png",pGRData]*)
(*Export["../graphs/GrowthRate3D.png",pDensGR]*)


(* ::Section:: *)
(*Precession frequency associated with the fastest mode *)


(* ::Subsubsection:: *)
(*Import data*)


(* ::Subsubsection:: *)
(*Precession frequency*)


(* ::Input:: *)
(*tabqx=Import["../data/Dump_Growth_Rate_Precession.hf5",{"Datasets","tabqx"}];*)
(*tabRot=Import["../data/Dump_Growth_Rate_Precession.hf5",{"Datasets","tabRot"}]; *)
(**)
(*RotTable=Partition[{tabqx,tabRot}//Transpose//Flatten,3];*)
(*a=Log10[RotTable[[All,1]]];*)
(*b=RotTable[[All,2]];*)
(*RotTable[[All,2]]=a;*)
(*RotTable[[All,1]]=b;*)
(**)


(* ::Subsubsection:: *)
(*Plot data*)


(* ::Input:: *)
(*(*Generic options for the plots.*)*)
(*fontname = "Helvetica";*)
(*fontsize = 15;*)
(*basestyleplot = {FontFamily -> fontname,FontSize->fontsize};*)
(*labelstyleplot = {FontFamily -> fontname, FontSize -> fontsize};*)
(*imagesize = Large;(*350;*)*)
(*(*Making the plot.*)*)
(*(*Range of the plot.*)*)
(*qmin = 0.01;*)
(*qmax =1;*)
(*xmin = 0.0;*)
(*xmax = 0.4;*)
(*zmin =0;*)
(*zmax = 5;*)
(*(*Axes of the plot.*)*)
(*ax = Row[{Style["q = ",fontname,fontsize],Style[FractionBox["Mdisk","Mdisk+MDH"]//DisplayForm,fontname,fontsize]}];*)
(*ay =  Row[{Style["x = ",fontname,fontsize],Style[FractionBox["Mbulb","Mdisk+Mbulb"]//DisplayForm,fontname,fontsize]}];*)
(*(*.*)*)
(*axstyle = {Black,Thickness[0.003]};*)
(*tickstyle = {Black,Thickness[0.003]};*)
(*(*.*)*)
(*qticks =  {{-3,"-3.0"},{-2.75,""},{-2.5,"-2.5"},{-2.25,""},{-2,"-2.0"},{-1.75,""},{-1.5,"-1.5"},{-1.25,""},{-1,"-1.0"},{-0.75,""},{-0.5,"-0.5"},{-0.25,""},{0,"0.0"}};*)
(*xticks =  {{0.0,"0.0"},{0.1,"0.1"},{0.2,"0.2"},{0.3,"0.3"},{0.4,"0.4"},{0.5,"0.5"},{0.6,"0.6"},{0.7,"0.7"},{0.8,"0.8"},{0.9,"0.9"},{1.0,"1.0"}};*)
(*zticks={1,2,3,4,5,6,7,8,9,10};RotTicks={1,2,3,4,5,6,7,8,9,10};*)


(* ::Input:: *)
(*pRot=ListContourPlot[RotTable,Contours->100,PlotRange->{{xmin,xmax},Log10[{qmin,qmax}],{zmin,zmax}},ColorFunction->"Rainbow",FrameLabel->{ay,ax},BaseStyle->basestyleplot,LabelStyle->labelstyleplot,PlotLegends->Placed[BarLegend[{(ColorData["Rainbow"][#]&),{zmin,zmax}},LegendLabel->"Rotation frequency",Method->{FrameStyle->axstyle,TicksStyle->tickstyle,Ticks->RotTicks}],Above],ImageSize->imagesize,FrameTicks->{qticks,xticks},FrameStyle->{axstyle,axstyle,axstyle,axstyle},PlotRangePadding->None]*)
(*pRot=Rasterize[pRot,RasterSize->1000];*)
(*pDensRot=ListPlot3D[RotTable,PlotRange->{{xmin,xmax},Log10[{qmin,qmax}],{zmin,zmax}},ColorFunction->"Rainbow",BaseStyle->basestyleplot,LabelStyle->labelstyleplot,ImageSize->imagesize,PlotRangePadding->None,AxesLabel->{ay,ax,"GR"},MeshFunctions->{#3&}];*)
(*pDensRot=Rasterize[pDensRot,RasterSize->1000]*)


(* ::Subsubsection:: *)
(*Save data*)


(* ::Input:: *)
(*Export["../graphs/PrecessionFrequency.png",pRot]*)
(*Export["../graphs/PrecessionFrequency3D.png",pDensRot]*)
