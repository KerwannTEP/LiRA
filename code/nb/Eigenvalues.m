(* ::Package:: *)

(* ::Section:: *)
(*Eigenvalues *)


(* ::Subsection:: *)
(*Import data*)


(* ::Input:: *)
(*tabEigVals=Import["../data/Dump_Eigenvalues.hf5",{"Datasets","tabEigVals"}];*)
(*tabEigValsPhys=Import["../data/Dump_Eigenvalues.hf5",{"Datasets","tabEigValsPhys"}];*)
(*tabEigValsReal=tabEigVals[[All,1]];*)
(*tabEigValsImag=tabEigVals[[All,2]];*)
(*tabEigValsPhysReal=tabEigValsPhys[[All,1]];*)
(*tabEigValsPhysImag=tabEigValsPhys[[All,2]];*)
(*tabEigValsPlane={tabEigValsReal,tabEigValsImag}//Transpose;*)
(*tabEigValsPhysPlane={tabEigValsPhysReal,tabEigValsPhysImag}//Transpose;*)
(**)


(* ::Subsection:: *)
(*Plot Nyquist diagramm*)


(* ::Input:: *)
(*imMax=2*)
(*pEig=ListPlot[{tabEigValsPlane},ImageSize->Large,PlotRange->All]*)
(*pEigPhys=ListPlot[{tabEigValsPhysPlane},ImageSize->Large,PlotRange->All]*)


(* ::Subsection:: *)
(*Save data*)


(* ::Input:: *)
(*Export["../graphs/Eigenvalues.png",pEig]*)
(*Export["../graphs/Physical_Eigenvalues.png",pEigPhys]*)
