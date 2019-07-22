(* ::Package:: *)

BeginPackage["OperatorGB`"]


Global`OperatorGB::usage="A few definitions need to be made by the user. Type ?SetUpRing for more information.";


CoeffQ::usage="CoeffQ[k_] should check whether k is in the base ring of the polynomial ring.";
Prod::usage="Prod[m1,m2,...] represents the non-commutative multiplication.";
ToProd::usage="Convert a polynomial from the built in non-commutative multiplication to the Prod data structure";
ToNonCommutativeMultiply::usage="Convert a polynomial from the Prod data structure to the built in non-commutative multiplication";


SetUpRing::usage="SetUpRing defines the non-commutative polynomial ring and a monomial order.

This method can be called with one list X = {x1,x2,...,xn} containing all indeterminates the user wants to work with. This then sets up the non-commutative polynomial ring \!\(\*TemplateBox[{},\n\"Rationals\"]\)<X> of non-commutative polynomials in the variables x1,x2,...,xn with rational numbers as coefficients. By default, this also sets up a graded lexicographic order where the indeterminates are ordered in increasing order according to their appearance in the list X, i.e. x1 < x2 < ... < xn.

The user has the option to change the monomial order. The second pre-defined order is a weighted graded lexicographic order. To switch to this order, the user has execute the command SortedQ := WeightedDegLex. Then, additionally, a list of weights {w1,w2,...} named Weight has to be provided. Type ?Weight for more information.

To define an individual order, the user can provide a binary function SortedQ[a,b] working on the set of all words that can be built from the alphabet of indeterminates in X. Type ?SortedQ for more information.

It is also possible to set up a multigraded lexicographic order. To do this, the user has to call SetUpRing with two lists X = {x1,x2,...,xn} and Y = {y1,y2,...,ym} as input. This defines the non-commutative polynomial ring \!\(\*TemplateBox[{},\n\"Rationals\"]\)<X,Y> of non-commutative polynomials in the variables x1,x2,...,xn,y1,y2,...,ym with rational numbers as coefficients. In this case, by default a multigraded lexicographic order is defined, where two monomials M1 and M2 are first compared by their degree only in indeterminates from Y, then by their degree only in indeterminates from X and finally, to break ties, a graded lexicographic order x1 < ... < xn < y1 < ... < ym is used. The multigraded lexicographic order is denoted by x1 < x2 < ... < xn << y1 < y2 < ... < ym.
"


Groebner::usage="Groebner[cofactors_,ideal_, maxiter:_?IntegerQ:10, OptionsPattern[{MaxDeg->Infinity,Info->False,Parallel->True,Sorted->False,OutputProd->False}]] executes at most maxiter iterations of the Buchberger algorithm to compute
a (partial) Groebner basis of an ideal. Additionally, for every new element in the Groebner basis a list
of cofactors is saved in the list cofactors forming a linear combination of the new element."


GroebnerWithoutCofactors::usage="GroebnerWithoutCofactors[ideal_,maxiter:_?IntegerQ:10,OptionsPattern[{MaxDeg->Infinity,Info->False,Parallel->True,Sorted->False,Criterion->False}]] executes at most maxiter iterations 
of the Buchberger algorithm to compute a (partial) Groebner basis of an ideal."


MultiplyOut::usage="To multiply out a list of cofactors given in terms of the built in non-commutative multiplication."


Rewrite::usage="Rewrite[vars, cofactors] rewrites a linear combination, which is safed in vars, with the elements from cofactors."


adj::usage="adj[A] represents the adjoint of the operator A"


ReducedForm::usage="ReducedForm[cofactors,G,exp] reduces the expression exp by the elements of G and saves the cofactors of the reduction process in the list cofactors. The argument exp can also
be a list of expressions, then all expressions are reduced."


ApplyRules::usage="ApplyRules[exp,G] reduces the expression exp using polynomials from the set G."


LeadingTerm::usage="Leading term of the polynomial w.r.t. the specified monomial ordering."
DegLex::usage="Degree Lexicographic order"
WeightedDegLex::usage="Weighted Degree Lexicographic order"
Weight::usage="A list {w1,w2,...} defining weights on the variables for the weighted graded lexicographic order. Then the weight w1 corresponds to the first variable in the input list of SetUpRing, w2 to the second one and so on."
SortedQ::usage="SortedQ[M1,M2] is a binary function working on the set of all words that can be built from the alphabet of variables given in SetUpRing. It decides, when given two words M1 and M2 as input, which of them is larger. SortedQ[M1,M2] returns True if M1 \[LessEqual] M2 and False otherwise.

M1 and M2 have to be given in form of lists containig only elements from the list(s) X (and Y) that where used in the call of the method SetUpRing. For example, if X = {x,y,z} then M1 could be of the form {x,x,y,z} and M2 could be {y,z,x,y}.
"


Quiver::usage="Data structure of a Quiver"
QSignature::usage="QSignature[poly,Q] returns the signature of the polynomial poly w.r.t. the quiver Q (not necessarily with unique lables)"
PlotQuiver::usage="Plot a quiver Q"


Certify::usage="Certifies whether a certain claim is a consequence of some assumptions via Groebner 
basis computations. Additionally, compatibility with a given quiver is checked."


(*Begin["`Private`"]*)


(* ::Subsection::Closed:: *)
(*Non-commutative multiplication*)


CoeffQ[_Rational]:=True
CoeffQ[_Integer]:=True
CoeffQ[_]:=False


Prod[a___,0,b___]=0;
Prod[a___,Prod[b___],c___]:=Prod[a,b,c]
Prod[a___,b_Plus,c___]:=(Prod[a,#,c]&/@b)
Prod[a___,c_?(CoeffQ[#]&),b___]:=c Prod[a,b]
Prod[a___,(d_?CoeffQ)*b_,c___]:=d Prod[a,b,c]
Prod[a___,b_List,c___]:=Prod[a,Sequence@@b,c]


ToProd[poly_]:= Expand[Prod[(poly//.NonCommutativeMultiply->Prod)]]


ToNonCommutativeMultiply[poly_]:= poly//.{Prod[]->1, Prod[a_]->Times[a], Prod[a_,b__]->NonCommutativeMultiply[a,b]}


(* ::Subsection::Closed:: *)
(*Setting up the ring*)


SetUpRing[vars_List]:= (
WordOrder = vars;
SortedQ := DegLex;
Print[Sequence@@Map[ToString[#,StandardForm]<>" < "&,vars[[;;-2]]] <> ToString[vars[[-1]],StandardForm]];
);


SetUpRing[knowns_List,unknowns_List]:= Module[{string},
	WordOrder = Join[knowns,unknowns]; 
	Knowns = knowns; 
	Unknowns = unknowns;
	SortedQ := MultiLex;
	
	string = Sequence@@Map[ToString[#,StandardForm]<>" < "&,knowns[[;;-2]]] <> ToString[knowns[[-1]],StandardForm];
	string = string <> " << ";
	string = string <> Sequence@@Map[ToString[#,StandardForm]<>" < "&,unknowns[[;;-2]]] <> ToString[unknowns[[-1]],StandardForm];
	Print[string];
]


(* ::Subsection::Closed:: *)
(*Ordering*)


(* ::Text:: *)
(*LeadingTerm[poly_] returns the leading term of poly (with regard to the term order set by the user). The return value is a triple {lc, lm, w}, where lc is the leading coefficient, lm is the leading monomial and w is the word corresponding to the leading monomial. *)


LeadingTerm[poly_Prod]:=
	List[1,List@@poly]


LeadingTerm[poly_]:=
	If[Head[Expand[poly]] === Times,
		List[First[poly],List@@Last[poly]],
		FindMaxTerm[Expand[poly]]
	]


FindMaxTerm[poly_]:= 
Module[{terms},
	terms = MonomialList[poly];
	Last[Sort[If[Head[#]===Times,
					List[#[[1]],List@@(#[[2]])],
					List[1,List@@#]]&/@terms,SortedQ[#1[[2]],#2[[2]]]&]]
]


(* ::Text:: *)
(*Multigraded Lexicographic order*)


MultiLex[a_List,b_List]:= 
Module[{V1a,V2a,V1b,V2b},
	V1a = Count[a,Alternatives@@Unknowns];
	V2a = Count[a,Alternatives@@Knowns];
	V1b = Count[b,Alternatives@@Unknowns];
	V2b = Count[b,Alternatives@@Knowns];

	If[(V1a < V1b) || (V1a === V1b && V2a < V2b), Return[True]];
	If[(V1a > V1b) || (V1a === V1b && V2a > V2b), Return[False]];
	DegLex[a,b]
]


(* ::Text:: *)
(*Degree Lexicographic order *)


DegLex[a_List,b_List]:=
Module[{i,l},
	l = Length[a];
	i = 1;
	If[l < Length[b],
		True, 
		If[l > Length[b],
			False,
			While[i <= l && Position[WordOrder,a[[i]]][[1,1]] === Position[WordOrder,b[[i]]][[1,1]],i++];
			If[i > l, True, Position[WordOrder,a[[i]]][[1,1]] < Position[WordOrder,b[[i]]][[1,1]]]
		]
	]
]


(* ::Text:: *)
(*Weighted Degree Lexicographic order. The weight has to be set by the user.*)


WeightedDegLex[a_List,b_List]:= 
Module[{i,la,lb,l},
	l = Length[a];
	la = WeightedLength[a];
	lb = WeightedLength[b];
	i = 1;
	If[la < lb,
		True,
		If[la > lb,
			False,
			While[i <= l && Position[WordOrder,a[[i]]][[1,1]] === Position[WordOrder,b[[i]]][[1,1]],i++];
			If[i > l, True, Position[WordOrder,a[[i]]][[1,1]] < Position[WordOrder,b[[i]]][[1,1]]]
		]
	]
]


WeightedLength[term_List]:=
Module[{i},
	Total[Weight[[#]]&/@ Table[Position[WordOrder,term[[i]]][[1,1]],{i,1,Length[term]}]]
]


(* ::Text:: *)
(*SortedQ[a_,b_] is a binary function working on the set of all words that can be built from the alphabet of modules defined in CyclicModules. It decides, when given two words as input, which of them is larger. *)
(*Hence, this function determines the monomial order. By defining an own SortedQ function the user can implement his own monomial order. By default the Degree Lexicographic order is set. It is also possible to switch to a weighted Degreelexicographic order by setting SortedQ := WeightedDegLex. Then the user also has to define a list Weight = {w1,w2,...} where w1 defines the weight for the first module in WordOrder, w2 defines the weight for the second module in WordOrder and so on.*)


SortedQ := DegLex;


(* ::Subsection::Closed:: *)
(*Ambiguities*)


(* ::Text:: *)
(*Reducible Words*)


ExtractReducibleWords[sys_]:=
	MapIndexed[{#1[[1]],#2[[1]]}&,sys]


(* ::Text:: *)
(*Overlaps*)


GenerateOverlaps[l1_List,l2_List,OptionsPattern[Parallel->True]]:=
	If[OptionValue[Parallel],
		Flatten[Parallelize[Outer[{Overlap[#1,#2],Overlap[#2,#1]}&,l1,l2,1],DistributedContexts->Automatic]],
		Flatten[Outer[{Overlap[#1,#2],Overlap[#2,#1]}&,l1,l2,1]]
	]


GenerateOverlaps[l_List,OptionsPattern[Parallel->True]]:= 
	If[OptionValue[Parallel] && Length[l] > 60,
		Flatten[Parallelize[Outer[Overlap,l,l,1],DistributedContexts->Automatic]],
		Flatten[Outer[Overlap,l,l,1]]
	]


Overlap[{v_List,i_Integer},{w_List,j_Integer}]:=
Module[{k,min},
	min = Min[Length[v],Length[w]];
	Reap[For[k=1,k<min,k++,
		If[Take[v,-k]===Take[w,k],
			Sow[Overlap[Join[v,Drop[w,k]],Drop[w,k],Drop[v,-k],{i,j}]]
		]];
	][[2]]
]


(* ::Text:: *)
(*Inclusions*)


GenerateInclusions[l_List,OptionsPattern[Parallel->True]]:= 
	If[OptionValue[Parallel] && Length[l] > 60,
		Flatten[Parallelize[Outer[Inclusion,l,l,1],DistributedContexts->Automatic]],
		Flatten[Outer[Inclusion,l,l,1]]
	]


GenerateInclusions[l1_List,l2_List,OptionsPattern[Parallel->True]]:= 
	If[OptionValue[Parallel],
		Flatten[Parallelize[Outer[{Inclusion[#1,#2],Inclusion[#2,#1]}&,l1,l2,1],DistributedContexts->Automatic]],
		Flatten[Outer[{Inclusion[#1,#2],Inclusion[#2,#1]}&,l1,l2,1]]
	]


Inclusion[{v_List,i_Integer},{w_List,j_Integer}]:=
Module[{k},
	Reap[If[Length[w]<Length[v],
		For[k=1,k+Length[w]-1<=Length[v],k++,
			If[v[[k;;(k+Length[w]-1)]]===w,
				Sow[Inclusion[v,v[[1;;(k-1)]],v[[(k+Length[w]);;]],{i,j}]]
			]
		]
	];
	][[2]]
]


(* ::Text:: *)
(*All ambiguities*)


GenerateAmbiguities[l_List,newpart_List,maxdeg_,OptionsPattern[Parallel->True]]:= 
	Select[Join[GenerateOverlaps[l,newpart,Parallel->OptionValue[Parallel]],
		GenerateInclusions[l,newpart,Parallel->OptionValue[Parallel]],
		GenerateAmbiguities[newpart,maxdeg,Parallel->OptionValue[Parallel]]], Length[#[[1]]] <= maxdeg &]


GenerateAmbiguities[l_List,maxdeg_,OptionsPattern[Parallel->True]] := 
	Select[Join[GenerateOverlaps[l,Parallel->OptionValue[Parallel]],
		GenerateInclusions[l,Parallel->OptionValue[Parallel]]], Length[#[[1]]] <= maxdeg &]


(* ::Text:: *)
(*Approach from the PhD Thesis*)


DeleteRedundantPhD[amb_List]:= 
Module[{Os,B={},pairs,toDelete,rules,i,j,wi,wj,k},
	Do[
		Os = Select[amb,#[[4,2]]===k&];
		pairs = Subsets[Os,{2}];
	
		toDelete = Map[(#/.{{Inclusion[_,wi_,_,{i_,_}],Inclusion[_,wj_,_,{j_,_}]}/; (i===j && wi!=wj) -> If[SortedQ[wj,wi],#[[1]],#[[2]]],
							{Inclusion[_,wi_,_,{i_,_}],Inclusion[_,wj_,_,{j_,_}]}/; (i=!=j) -> If[i > j,#[[1]],#[[2]]],
							{Overlap[_,wi_,_,{i_,_}],Inclusion[_,wj_,_,{j_,_}]}/; (i===j && wi !=wj) -> If[SortedQ[wj,wi],#[[1]],#[[2]]],
							{Overlap[_,wi_,_,{i_,_}],Inclusion[_,wj_,_,{j_,_}]}/; (i=!=j) -> If[i > j,#[[1]],#[[2]]],
							
							#->Nothing}
		
		)&,pairs];
			
		Os = DeleteCases[Os,Alternatives@@toDelete];
		If[Length[toDelete] > 0,
			Print["Removing ", Length[toDelete], " ambiguities..."]
		];
		B = Join[B,Os];
	,{k,Max[amb[[All,4,2]]]}];
	B
]


(* ::Text:: *)
(*Approach from Mora*)


DeleteRedundant[ambInput_List,OptionsPattern[{Info->False}]]:= 
Module[{selected,i,amb,result,f,t},
	t = AbsoluteTiming[
	amb = SortBy[ambInput,Length[#[[1]]&]];
	result = {};
	Do[
		selected = Select[amb,Max[#[[4]]]===i &];
		While[Length[selected] > 0,
			f = First[selected];
			selected = Drop[selected,1];
			AppendTo[result,f];
			selected = DeleteCases[selected,_[{___,Sequence@@f[[1]],___},__]];
		];
		,{i,Min[Flatten[amb[[All,4]]]],Max[Flatten[amb[[All,4]]]]}
	];
	][[1]];
	If[OptionValue[Info],
		Print["Removed ", Length[amb] - Length[result], " ambiguities in ",t]];
	result
]


(* ::Subsection::Closed:: *)
(*S-Polynomials*)


(* ::Text:: *)
(*SPoly[amb,fi,fj] computes the S-polynomial corresponding to the ambiguity amb, which comes from the two reduction rules fi and fj. Additionally, the linear combination how the S-polynomial was computed from fi and fj is returned in a list. *)


SPoly[amb:_Overlap|_Inclusion,fi_,fj_]:=
Module[{A,C},
	C = Prod@@amb[[2]];
	A = Prod@@amb[[3]];
	If[amb[[0]]=== Overlap,
			(*Overlap[ABC,C,A]*)
			{Prod[fi[[2]],C] - Prod[A,fj[[2]]],
				{{A,ToPoly[fj],Prod[]},{-Prod[],ToPoly[fi],C}}},
			(*Inclusion[CBA,C,A]*)
			{fi[[2]] - Prod[C,fj[[2]],A],
			{{C,ToPoly[fj],A},{-Prod[],ToPoly[fi],Prod[]}}}
	]
]


(* ::Text:: *)
(*Same as Poly but without returning the linear combination.*)


SPoly2[amb:_Overlap|_Inclusion,fi_,fj_]:=
	If[amb[[0]]=== Overlap,
			(*Overlap[ABC,C,A]*)
			Prod[fi[[2]],amb[[2]]] - Prod[amb[[3]],fj[[2]]],
			(*Inclusion[CBA,C,A]*)
			fi[[2]] - Prod[amb[[2]],fj[[2]],amb[[3]]]
	]


(* ::Subsection::Closed:: *)
(*Gr\[ODoubleDot]bner basis*)


(* ::Text:: *)
(*Implementation of the Buchberger algorithm to compute a (partial) Groebner basis of the reduction system sys with at most maxiter iterations being executed (default: 10). The return value is a reduction system {f1,...,fn,g1,...gm} consisting of the elements f1,....,fn from sys and new elements g1,...,gm. For every new element g, a pair {g, l} is saved in the list cofactors, where l is a list forming a linear combination of g consisting of elements from sys and certain cofactors.*)
(**)
(*OptionPattern:*)
(*	- Ignore (default: 0): A non-negative integer that determines how many elements of the input will be ignored during the first computation of the ambiguities. *)
(*	- MaxDeg (default: Infinity): Only ambiguities with degree smaller than or equal to MaxDeg will be considered during the Groebner basis computation (larger ambiguities are simply ignored). *)
(*	- Info (default: False): Prints information about the computation progress.*)
(*	- Parallel (default: True): Determines whether the computations for which it is possible, are executed in parallel (which speeds up the computation) or in series.*)
(*	- Sorted (default: False):  Sorts the ambiguities before processing in ascending order. This speeds up the computation but results in a different (partial) Groebner basis.*)
(*	- OutputProd (default: False): If this OptionPattern is set to True, the output, i.e. the Groebner basis and the list of cofactors, is given in the Prod data structure. Otherwise, Mathematica's*)
(*	non-commutative multiplication is used.*)
(**)


SetAttributes[Groebner,HoldFirst]

Groebner[cofactors_,ideal_, maxiter:_?IntegerQ:10, OptionsPattern[{Criterion->False,Ignore->0,MaxDeg->Infinity,Info->False,Parallel->True,Sorted->False,OutputProd->False,Rewrite->True}]]:=
Module[{count,spol,lt,info,p,h,G,r,t1,t2,rules,sorted,oldlength,parallel,hrule,maxdeg,outputProd,criterion,i},
	info = OptionValue[Info];
	sorted = OptionValue[Sorted];
	parallel = OptionValue[Parallel];
	maxdeg = OptionValue[MaxDeg];
	outputProd = OptionValue[OutputProd];
	criterion = OptionValue[Criterion];

	If[Head[cofactors]=!=List,cofactors={}];
	
	G = CreateRedSys[ideal];
	oldlength = Length[G];
	t1 = 0; t2 = 0; count = 0;
	If[info,Print["G has ", Length[G]," elements in the beginning."],Print[]];

	spol = CheckResolvability[G,OptionValue[Ignore],Criterion->criterion,MaxDeg->maxdeg,Info->info,Sorted->sorted,Parallel->parallel];
	rules = ExtractRules[G];

	While[Length[spol] > 0 && count < maxiter,
		t1 = AbsoluteTiming[
		i = Length[spol];
		Monitor[Do[
			(*reduce it*)
			r = Reap[p[[1]]//.rules]; 
			h = r[[1]];
			If[h =!= 0,
				If[Length[r[[2]]] > 0,
					p[[2]]\[NonBreakingSpace]= Join[p[[2]],r[[2,1]]]
				];
				lt = LeadingTerm[h];
				If[lt[[1]] =!= 1, 
					h = Expand[1/lt[[1]]*h]; p[[2]] = (ReplacePart[#,1 -> 1/lt[[1]]*#[[1]]]&/@ p[[2]])
				];
				hrule = CreateRedSys[h];
				AppendTo[G,hrule];
				AppendTo[cofactors,{h,p[[2]]}]; 
				AppendTo[rules,Sequence@@ExtractRules[{hrule}]];
			];
			i--;
		,{p,spol}];,i];][[1]];

		If[info, Print["The second reduction took ", t1]];
		count++;
		If[info,Print["Iteration ",count, " finished. G has now ", Length[G]," elements\n"]];
		If[count < maxiter, 
			spol = CheckResolvability[G,oldlength,Criterion->criterion,MaxDeg->maxdeg,Info->info,Sorted->sorted,Parallel->parallel];
			oldlength = Length[G];
		];
	];
	
	G = ToPoly[G];
	
	If[OptionValue[Rewrite],
		If[info, Print["Rewriting the cofactors has started."]];
		t2 = AbsoluteTiming[
			RewriteGroebner[cofactors,Info->info,OutputProd->outputProd];
		][[1]];
		If[info, Print["Rewriting the cofactors took in total ", t2]];
	];
	If[outputProd,
		G,
		ToNonCommutativeMultiply[G]
	]
]


(* ::Text:: *)
(*CheckResolvability[sys,Info->False,Parallel->True,Sorted->False] returns all S-polynomials from the reduction system sys which can not be reduced to zero. Additionally, *)
(*for each S-polynomial a list containing the linear combination how the S-polynomial was generated from the elements of sys is returned.*)
(*For a description of the OptionPatterns see the documentation of the Groebner method.*)


CheckResolvability[sys_,oldlength:_?IntegerQ:0,OptionsPattern[{Criterion->False,MaxDeg->Infinity,Info->False,Parallel->True,Sorted->False}]]:=
Module[{amb,spol,info,rules,parallel,words,r,t1,t2,t3},
	info = OptionValue[Info];
	parallel = OptionValue[Parallel];

	(*generate ambiguities*)
	words = ExtractReducibleWords[sys];
	t1 = AbsoluteTiming[
		amb = GenerateAmbiguities[words[[;;oldlength]],words[[oldlength+1;;]],OptionValue[MaxDeg],Parallel->parallel];
	][[1]];
	If[info,Print[Length[amb]," ambiguities in total (computation took ",t1, ")"]];
	
	(*process ambiguities*)
	If[OptionValue[Criterion],
		amb = DeleteRedundant[amb,Info->info]
	];
	If[OptionValue[Sorted],amb = Sort[amb]];
	
	(*generate S-polynomials*)
	t2 = AbsoluteTiming[
	spol = DeleteCases[
		If[parallel,
			ParallelMap[SPoly[#,sys[[#[[4,1]]]],sys[[#[[4,2]]]]]&,amb,DistributedContexts->Automatic,Method->"CoarsestGrained"],
			Map[SPoly[#,sys[[#[[4,1]]]],sys[[#[[4,2]]]]]&,amb]
		]
	,{0,___}];
	][[1]];
	If[info,Print["Generating S-polys: ",t2 ," (",Length[spol], " in total)"]];
	
	(*reduce S-polynomials*)
	rules = ExtractRules[sys];
	(*parallelizing this makes it only slower*)
	t3 = AbsoluteTiming[
	spol = DeleteDuplicatesBy[DeleteCases[Map[(r = Reap[#[[1]]//.rules]; If[r[[1]]=!= 0,
											If[Length[r[[2]]] > 0,
												{r[[1]],Join[#[[2]],r[[2,1]]]},
												{r[[1]],#[[2]]}
												],
											{}
											])&,spol],{}],First];
	][[1]];
	If[info,Print["Reducing S-polys: ",t3, " (",Length[spol], " remaining)"]];
	spol
]


(* ::Text:: *)
(*ExtractRules[vars,sys] generates a set of reduction rules out of the reduction system sys with the special feature that the cofactors of the reduction will be sowed whenever such a rule is applied. They can then be reaped using the Reap command.*)
(*Sys has to consist of pairs {word,func}, which can be obtained using the method CreateRedSys.*)


ExtractRules[sys_]:=
Module[{a,b,coeff,i,p,q},
	a=Unique[];b=Unique[];coeff=Unique[];
	Table[
		p = Prod[Pattern[Evaluate[a],BlankNullSequence[]],i[[1]],Pattern[Evaluate[b],BlankNullSequence[]]];
		q = Expand[Evaluate[coeff]*Prod[Evaluate[a],i[[2]],Evaluate[b]]];
		With[{x ={-Evaluate[coeff]*Prod[Evaluate[a]],ToPoly[i],Prod[Evaluate[b]]},y=q},Alternatives[Pattern[Evaluate[coeff],BlankNullSequence[]]*p,p]:> Hold[Sow[x];y]]
	,{i,sys}]//ReleaseHold
]


SetAttributes[RewriteGroebner,HoldFirst]

RewriteGroebner[cofactors_,OptionsPattern[{Info->False,OutputProd->False}]]:=
Module[{a,b,i,j,count,rules,info,occurring},
	info = OptionValue[Info];
	rules = {};
	count = Length[cofactors];
	
	Monitor[Do[
		(*only choose those rules which will actually be used*)
		occurring = Position[cofactors[[;;i,1]],Alternatives@@DeleteDuplicates[cofactors[[i+1,2,All,2]]],{1}];
		rules = Map[{a__,#[[1]],b__} -> Sequence@@Table[{Prod[Evaluate[a],j[[1]]],j[[2]],Prod[j[[3]],Evaluate[b]]},{j,#[[2]]}]&,Extract[cofactors,occurring]];
		(*reduce the cofactorlist*)
		cofactors[[i+1,2]] = cofactors[[i+1,2]]/.rules//CollectLeft//ExpandLeft;
		count = count - 1;
	,{i,Length[cofactors]-1}];
	,count];
	If[!OptionValue[OutputProd],
		cofactors = ToNonCommutativeMultiply[cofactors];
	];
]


(* ::Subsection::Closed:: *)
(*Cofactor criterion*)


(* ::Text:: *)
(*Does f divide p?*)


Divides[p_,f_,i_]:= 
Module[{termsP,termsF,x,y},
	termsP = Map[If[!CoeffQ[#[[1]]],{1,#},#]&,(MonomialList[p]/.Times->List)/.Prod->List];
	termsF = Map[If[!CoeffQ[#[[1]]],{1,#},#]&,(MonomialList[f]/.Times->List)/.Prod->List];
	If[MatchQ[termsP,Map[{#[[1]],{x___,Sequence@@#[[2]],y___}}&,termsF]],Print[i];{i},Nothing]
]


CofactorCriterion[spol_,cofactors_]:=
Module[{cofactorPolies,toDelete,toDelete1, toDelete2, spolTerms,t,t1,t2,spolPolies},
	t = AbsoluteTiming[
	spolPolies = spol;
	cofactorPolies = DeleteDuplicates[Flatten[Map[#/.List->Prod&,cofactors,{2}]]];
	][[1]];
	Print["Cofactors to Polies took ", t];	
	t1 = AbsoluteTiming[
	MakeMonic[spolPolies];
	MakeMonic[cofactorPolies];
	Print["All monic"];
	toDelete1 = Outer[Divides,spolPolies,cofactorPolies,Range[Length[spol]]];
	Print["Length[toDelete1] = ", Length[toDelete1]];
	][[1]];
	Print["Finding them took ", t1];
	toDelete1
]


(* ::Subsection::Closed:: *)
(*F4*)


SetAttributes[F4,HoldFirst];

F4[cofactors_,ideal_, maxiter:_?IntegerQ:10, OptionsPattern[{N->50,Criterion->False,Ignore->0,MaxDeg->Infinity,Info->False,Parallel->True,Sorted->False,OutputProd->False,Rewrite->True}]]:=
Module[{count,spol,lt,info,G,t1,t2,sorted,oldlength,parallel,maxdeg,n,L,cofactorsL,lc,a,b,rules},
	info = OptionValue[Info];
	sorted = OptionValue[Sorted];
	parallel = OptionValue[Parallel];
	maxdeg = OptionValue[MaxDeg];
	
	cofactors = {};

	G = ToProd/@ideal;
	lc = MakeMonic[G];

	oldlength = Length[G];
	t1 = 0; t2 = 0; count = 0;
	If[info,Print["G has ", Length[G]," elements in the beginning."],Print[]];

	spol = CheckResolvabilityF4[G,OptionValue[Ignore],Criterion->OptionValue[Criterion],MaxDeg->maxdeg,Info->info,Sorted->sorted,Parallel->parallel];
	
	While[count < maxiter,
		t1 = AbsoluteTiming[
		Monitor[While[Length[spol]\[NonBreakingSpace]> 0,
			(*pick some S-polynomials and reduce them*)
			n = Min[OptionValue[N],Length[spol]];
			{L,cofactorsL}=Reduction[spol[[;;n]],G];
			spol = Drop[spol,n];
			G = Join[G,L];
			cofactors = Join[cofactors,MapThread[{#1,#2}&,{L,cofactorsL}]];
		];
		,Length[spol]];][[1]];
		
		If[Length[G]===oldlength,Print["All S-polynomials could be reduced to 0."];Break[]];

		If[info, Print["The reduction took ", t1]];
		count++;
		If[info,Print["Iteration ",count, " finished. G has now ", Length[G]," elements\n"]];
		If[count < maxiter, 
			spol = CheckResolvabilityF4[G,oldlength,Criterion->OptionValue[Criterion],MaxDeg->maxdeg,Info->info,Sorted->sorted,Parallel->parallel];
			oldlength = Length[G];
		];
	];
	
	If[info, Print["Rewriting the cofactors has started."]];
	t2 = AbsoluteTiming[
	(*take care of leading coefficients*)
	rules = MapIndexed[{a_,#1,b_}->{a/lc[[First[#2]]],Expand[lc[[First[#2]]]*#1],b}&,G[[;;Length[ideal]]]];
	cofactors = cofactors/.rules;
	RewriteGroebner[cofactors,Info->info];
	][[1]];
	If[info, Print["Rewriting the cofactors took in total ", t2]];

	ToNonCommutativeMultiply[G]
]


(* ::Subsubsection:: *)
(*Symbolic Preprocessing & Reduction*)


SetAttributes[SymbolicPreprocessing,HoldFirst]

SymbolicPreprocessing[cofactorsF_,L_,G_]:=
Module[{F,T,lt,g,rules,a,b},
	F = L;
	T = DeleteDuplicates[Flatten[Monomials/@F,1]];
	lt = (LeadingTerm/@G)[[All,2]];
	rules = Table[With[{x = i},{a___,Sequence@@lt[[i]],b___}:>(AppendTo[cofactorsF,{Prod[a],G[[x]],Prod[b]}];Prod[a,G[[x]],b])],{i,Length[G]}];
	While[Length[T] > 0,
		g = DeleteCases[T/.rules,_List];
		F = Join[F,g];
		T = Complement[DeleteDuplicates[Flatten[Monomials/@g,1]],T];
	];
	lt = MakeMonic[F];
	cofactorsF = MapIndexed[{#1[[1]]/lt[[#2[[1]]]],#1[[2]],#1[[3]]}&,cofactorsF];
	cofactorsF = DeleteDuplicatesBy[cofactorsF,Prod[#]&];
	DeleteDuplicates[F]
]


Reduction[L_,G_]:=
Module[{F,M,lt,columns,FPlus,a,ct1,t2,t3,t4,cofactorsF,A,cofactors,pos},
	t1 = AbsoluteTiming[
	cofactorsF = L[[All,2]];
	F = SymbolicPreprocessing[cofactorsF,L[[All,1]],G];
	][[1]];
	
	t2 = AbsoluteTiming[
	lt = (LeadingTerm/@F)[[All,2]];
	(*sort in descending order*)
	columns = Reverse[Sort[DeleteDuplicates[Flatten[Monomials/@F,1]],SortedQ]];
	M = SparseArray[Flatten[MapIndexed[#1/.{Plus->List,c_*Prod[a___]->({#2[[1]],{a}}->c),Prod[a___]->({#2[[1]],{a}}->1)}&,F]/.MapIndexed[#1->#2[[1]]&,columns]]];
	][[1]];
	
	t3 = AbsoluteTiming[
	{A,M} = HermiteDecomposition[M];
	FPlus = DeleteCases[M.(ToProd/@columns),0];
	][[1]];
	
	t4 = AbsoluteTiming[
	cofactors = Map[(Table[{#[[i]]*cofactorsF[[i,1]],cofactorsF[[i,2]],cofactorsF[[i,3]]},{i,Length[F]}]//ExpandLeft)&,A];
	cofactors = DeleteCases[cofactors,{}];
	pos = Position[FPlus,p_/;!MemberQ[lt,LeadingTerm[p][[2]]],{1},Heads->False]; 
	][[1]];
	Print["SymPre: ", t1,", Setup: ", t2, ", RRed: ", t3, ", Cofactor stuff: ", t4];
	{Extract[FPlus,pos],Extract[cofactors,pos]}
]


(* ::Subsubsection:: *)
(*Additional stuff*)


Monomials[f_]:= (MonomialList[f]/.Times[_,a_]->a)/.Prod->List


(* ::Text:: *)
(*Input format: normal polynomials but in Prod data structure*)


SPolyF4[amb:_Overlap|_Inclusion,fi_,fj_]:=
	If[amb[[0]]=== Overlap,
			(*Overlap[ABC,C,A]*)
			Sequence@@{{Prod[fi,amb[[2]]],{Prod[],fi,Prod[amb[[2]]]}},{Prod[amb[[3]],fj],{Prod[amb[[3]]],fj,Prod[]}}},
			(*Inclusion[CBA,C,A]*)
			Sequence@@{{fi,{Prod[],fi,Prod[]}},{Prod[amb[[2]],fj,amb[[3]]],{Prod[amb[[2]]],fj,Prod[amb[[3]]]}}}
	]


CheckResolvabilityF4[G_,oldlength:_?IntegerQ:0,OptionsPattern[{Criterion->False,MaxDeg->Infinity,Info->False,Parallel->True,Sorted->False}]]:=
Module[{amb,spol,info,t1,t2,lt,sorted,maxdeg,parallel,words},

	info = OptionValue[Info];
	sorted = OptionValue[Sorted];
	parallel = OptionValue[Parallel];
	maxdeg = OptionValue[MaxDeg];

	(*generate ambiguities*)
	words = MapIndexed[{LeadingTerm[#1][[2]],First[#2]}&,G];
	t1 = AbsoluteTiming[
		amb = GenerateAmbiguities[words[[;;oldlength]],words[[oldlength+1;;]],maxdeg,Parallel->parallel];
	][[1]];
	If[info,Print[Length[amb]," ambiguities in total (computation took ",t1, ")"]];
	If[OptionValue[Criterion],
		amb = DeleteRedundant[amb,Info->info]
	];
	If[sorted,amb = Sort[amb]];
	
	(*generate S-polynomials*)
	t2 = AbsoluteTiming[
		If[parallel,
			spol = ParallelMap[SPolyF4[#,G[[#[[4,1]]]],G[[#[[4,2]]]]]&,amb,DistributedContexts->Automatic,Method->"CoarsestGrained"],
			spol = Map[SPolyF4[#,G[[#[[4,1]]]],G[[#[[4,2]]]]]&,amb]
		];		
	][[1]];
	spol = DeleteDuplicatesBy[spol,First];
	If[info,Print["Generating S-polys: ",t2]];
	If[info, Print[Length[spol]," critical polynomials were generated."]];
	spol
]


(* ::Subsection::Closed:: *)
(*Find Cofactors*)


(* ::Text:: *)
(*ReducedForm[cofactors,G,exp] can be used to reduce the expression exp with the elements of G. The linear combination of these reduction steps is saved in the list cofactors. *)
(*Exp can also be a list of expressions. *)


SetAttributes[ReducedForm,HoldFirst]

ReducedForm[cofactors_,G_,exp_]:=
Module[{t,rules},
	cofactors = {};
	rules = ExtractRules[CreateRedSys[G]];
	t = Reap[ToProd[exp]//.rules];
	If[Length[t[[2]]]>0,
		cofactors = ReplacePart[#,1->-#[[1]]]&/@t[[2,1]];
		cofactors = ToNonCommutativeMultiply[cofactors]
	];
	ToNonCommutativeMultiply[t[[1]]]
]


ReducedForm[cofactors_,G_,exp:_?ListQ]:=
Module[{i,t,rules},
	cofactors = Table[{},Length[exp]];
	rules = ExtractRules[CreateRedSys[G]];
	Table[
		t = Reap[ToProd[exp[[i]]]//.rules];
		If[Length[t[[2]]]>0,
			cofactors[[i]] = ReplacePart[#,1->-#[[1]]]&/@t[[2,1]];
			cofactors[[i]] = ToNonCommutativeMultiply[cofactors[[i]]];
		];
		ToNonCommutativeMultiply[t[[1]]]
	,{i,Length[exp]}
	]
]


(* ::Subsection::Closed:: *)
(*Gr\[ODoubleDot]bner basis without cofactors*)


(* ::Text:: *)
(*Implementation of the Buchberger algorithm to compute a (partial) Groebner basis of the reduction system sys with at most maxiter iterations being executed (default: 10). The return value is a reduction system {f1,...,fn,g1,...gm} consisting of the elements f1,....,fn from sys and new elements g1,...,gm. *)
(**)
(*OptionPattern:*)
(*	- MaxDeg (default: Infinity): Only ambiguities with degree smaller than or equal to MaxDeg will be considered during the Groebner basis computation (larger ambiguities are simply ignored). *)
(*	- Info (default: False): Prints information about the computation progress.*)
(*	- Parallel (default: True): Determines whether the computations for which it is possible, are executed in parallel (which speeds up the computation) or in series.*)
(*	- Sorted (default: False):  Sorts the ambiguities before processing in ascending order. This speeds up the computation but results in a different (partial) Groebner basis.*)
(*	-Criterion (default: False): Deletes potentially redundant elements of the Groebner basis during the computation. More specifically, if h | g_i for some g_i in G and spol(h,g_i = 0), then G = G\{g_i} u {h}. This criterion speeds up the computation but yields a different, in general smaller, (partial) Groebner basis.*)
(*	*)
(*	*)


GroebnerWithoutCofactors[ideal_,maxiter:_?IntegerQ:10,OptionsPattern[{Ignore->0, MaxDeg->Infinity,Info->False,Parallel->True,Sorted->False,Criterion->False}]]:=
Module[{count,spol,p,h,G,lt,info,t,rules,criterion,oldlength,maxdeg,incl,pos,sorted,parallel,syslength,i},

	info = OptionValue[Info];
	criterion = OptionValue[Criterion];
	sorted = OptionValue[Sorted];
	parallel = OptionValue[Parallel];
	maxdeg = OptionValue[MaxDeg];

	G = CreateRedSys[ideal];
	syslength = oldlength = Length[G];
	If[info,Print["G has ", Length[G]," elements in the beginning."];Print[]];
	count = 0; t = 0;

	spol = CheckResolvability2[G,OptionValue[Ignore],Criterion->criterion,MaxDeg->maxdeg,Info->info,Sorted->sorted,Parallel->parallel];
	rules = ExtractRules2[G];

	While[Length[spol] > 0 && count < maxiter,
		i = Length[spol];
		t = AbsoluteTiming[Monitor[Do[
			h = p//.rules;
			If[h =!= 0,  
				lt = LeadingTerm[h];
				If[lt[[1]] =!= 1, h = Expand[1/lt[[1]]*h]];
				AppendTo[G,CreateRedSys[h]];
				AppendTo[rules,Sequence@@ExtractRules2[{G[[-1]]}]]
			];
			i--;
		,{p,spol}];,i];][[1]];
		If[info, Print["The reduction took ", t]];
		count = count + 1;
		If[info,Print["Iteration ",count, " finished. G has now ", Length[G]," elements"];Print[]];
		If[count < maxiter,
			spol = CheckResolvability2[G,oldlength,Criterion->criterion,MaxDeg->maxdeg,Info->info,Sorted->sorted,Parallel->parallel];
			oldlength = Length[G];
		];
	];
	ToNonCommutativeMultiply[ToPoly[G]]
]


CheckResolvability2[sys_,oldlength:_?IntegerQ:0,OptionsPattern[{Criterion->False,MaxDeg->Infinity,Info->False,Sorted->False,Parallel->True}]]:=
Module[{amb,spol,info,t1,t2,lists,rules,words,parallel},
	info = OptionValue[Info];
	parallel = OptionValue[Parallel];

	words = ExtractReducibleWords[sys];
	rules = ExtractRules2[sys];
	
	(*generate ambiguities*)
	t1 = AbsoluteTiming[
		amb = GenerateAmbiguities[words[[;;oldlength]],words[[oldlength+1;;]],OptionValue[MaxDeg],Parallel->True]
	][[1]];
	If[info,Print[Length[amb]," ambiguities in total (computation took ", t1, ")"]];
	
	(*process ambiguities*)
	If[OptionValue[Criterion],
		amb = DeleteRedundant[amb,Info->info]
	];
	If[OptionValue[Sorted],amb = Sort[amb]];
	
	(*generate and reduce S-polynomials*)
	t2 = AbsoluteTiming[
	If[parallel && Length[amb] > 300,
		spol = DeleteDuplicates[DeleteCases[ParallelMap[SPoly2[#,sys[[#[[4,1]]]],sys[[#[[4,2]]]]]//.rules &,amb,DistributedContexts->Automatic,Method->"ItemsPerEvaluation" -> 1000],0]],
		spol = DeleteDuplicates[DeleteCases[Map[SPoly2[#,sys[[#[[4,1]]]],sys[[#[[4,2]]]]]//.rules &,amb],0]]];
	][[1]];

	If[info, Print[Length[spol]," different S-polynomials did not reduce to 0 (computation took ",t2,")"]];
	spol
]


(* ::Text:: *)
(*Parallelising this function does not make sense. It is faster when executed sequential.*)


ExtractRules2[sys_]:=
Module[{a,b,i},
	a=Unique[];b=Unique[];
	Table[Prod[Pattern[Evaluate[a],BlankNullSequence[]],i[[1]],Pattern[Evaluate[b],BlankNullSequence[]]]->
		Prod[a,i[[2]],b],{i,sys}]
]


ApplyRules[expr_List,G_]:= ApplyRules[#,G]&/@expr


ApplyRules[expr_,G_]:=
	ToNonCommutativeMultiply[ToProd[expr]//.ExtractRules[CreateRedSys[G]]]


(* ::Subsection::Closed:: *)
(*Additional stuff*)


NormalizePoly[poly_]:= 
	Expand[1/LeadingTerm[poly][[1]]*poly]


Remainder[poly_, lt:{_,_List}]:=
	poly - lt[[1]]*Prod@@lt[[2]]


CreateRedSys[polies_List]:=
	Map[CreateRedSys,polies]


CreateRedSys[p_]:=
Module[{lt,poly},
	poly = ToProd[p];
	lt = LeadingTerm[poly];
	{lt[[2]], Expand[-1/lt[[1]]*Remainder[poly,lt]]}
]


ToPoly[poly:{_List,_Prod|_Plus|_Times}]:= Prod@@poly[[1]] - poly[[2]]


ToPoly[sys_]:= ToPoly/@sys;


Rewrite[spolfactors:List[RepeatedNull[List[RepeatedNull[{__,__,__}]]]],cofactors_,OptionsPattern[InputProd->False]]:= 
	Map[Rewrite[#,cofactors,InputProd->OptionValue[InputProd]]&,spolfactors]


Rewrite[spolfactor_, cofactor_,OptionsPattern[InputProd->False]]:=
Module[{a,b,i,j,rules,occurring,result,spolfactors,cofactors},
	spolfactors = Map[ToProd,spolfactor,{2}];
	If[OptionValue[InputProd],
		cofactors = cofactor,
		cofactors = Map[{ToProd[#[[1]]],Map[ToProd,#[[2]],{2}]}&,cofactor];
	];
	If[cofactors === {},
		result = spolfactors,
		occurring = Position[cofactors[[All,1]],Alternatives@@DeleteDuplicates[spolfactors[[All,2]]],{1}];
		rules = Map[{a__,#[[1]],b__} -> Sequence@@Table[{Prod[a,#[[2,j,1]]],#[[2,j,2]],Prod[#[[2,j,3]],b]},{j,Length[#[[2]]]}]&,Extract[cofactors,occurring]];
		result = spolfactors/.rules;
	];
	If[OptionValue[InputProd],
		result//CollectLeft//ExpandLeft,
		ToNonCommutativeMultiply[result//CollectLeft//ExpandLeft]
	]
]


MultiplyOut[cofactors_List]:=Expand[ToNonCommutativeMultiply[Total[Map[ToProd,cofactors]]]]


SetAttributes[MakeMonic,HoldFirst];

MakeMonic[ideal_]:=
Module[{lc},
	lc = (LeadingTerm/@ideal)[[All,1]];
	ideal = Expand[ideal/lc];
	lc
]


Interreduce[ideal_,OptionsPattern[{InputProd->False}]]:=
Module[{G,rules,i,s,gi,cofactors,r,lt,sys,a,b,coeff,p,q,newPart},
	(*set everything up*)
	G = If[OptionValue[InputProd],
		ideal,
		ToProd/@ideal
	];
	sys = CreateRedSys[ideal];
	a=Unique[];b=Unique[];
	coeff = Unique[];
	rules = ReleaseHold[Table[
		p = Prod[Pattern[Evaluate[a],BlankNullSequence[]],sys[[i,1]],Pattern[Evaluate[b],BlankNullSequence[]]];
		q = Expand[Evaluate[coeff]*Prod[Evaluate[a],sys[[i,2]],Evaluate[b]]];
		With[{x ={-Evaluate[coeff]*Prod[Evaluate[a]],i,Prod[Evaluate[b]]},y=q},Alternatives[Pattern[Evaluate[coeff],BlankNullSequence[]]*p,p]:> Hold[Sow[x];y]]
	,{i,Length[sys]}]];
	i = 1; s = Length[G];
	cofactors = Map[{{Prod[],#,Prod[]}}&,G];
	
	(*actual interreduction*)
	While[i < s,
		If[G[[i]]===Null,i++;Continue[]];
		r = Reap[gi = G[[i]]//.Drop[rules,{i}]];
		If[gi===0,
			rules[[i]] = Null->Null;
			G[[i]] = Null;
			cofactors[[i]] = Null;
			i++,
			If[gi =!= G[[i]],
				lt = LeadingTerm[gi];
				r = r[[2,1]];
				newPart = Flatten[Table[Map[{Prod[r[[i,1]],#[[1]]],#[[2]],Prod[#[[3]],r[[i,3]]]}&,cofactors[[r[[i,2]]]]],{i,Length[r]}],1];
				cofactors[[i]] = Join[cofactors[[i]],newPart];
				If[lt[[1]] =!= 1, 
					gi = Expand[1/lt[[1]]*gi]; 
					cofactors[[i]] = (ReplacePart[#,1 -> 1/lt[[1]]*#[[1]]]&/@ cofactors[[i]]);
				];
				G[[i]]=gi;
				sys = CreateRedSys[gi];
				p = Prod[Pattern[Evaluate[a],BlankNullSequence[]],sys[[1]],Pattern[Evaluate[b],BlankNullSequence[]]];
				q = Expand[Evaluate[coeff]*Prod[Evaluate[a],sys[[2]],Evaluate[b]]];
				rules[[i]] = ReleaseHold[With[{x ={-Evaluate[coeff]*Prod[Evaluate[a]],i,Prod[Evaluate[b]]},y=q},Alternatives[Pattern[Evaluate[coeff],BlankNullSequence[]]*p,p]:> Hold[Sow[x];y]]];
				i = 1,
				i++;
			];
		];
	];
	{DeleteCases[G,Null],Map[ExpandLeft[CollectLeft[#]]&,DeleteCases[cofactors,Null]]}
]


(* ::Subsection::Closed:: *)
(*Collect cofactors*)


CollectLeft[triples_List]:=
Module[{x,y},
	Table[Prepend[x,Total[Cases[triples,{y_,Sequence@@x}:>y]]],{x,DeleteDuplicates[Rest/@triples]}]
]


ExpandLeft[triples_List]:=
Module[{x},
	Function[x,Sequence@@Which[
		x[[1]]===0,{},
		x[[1,0]]===Plus,Map[Prepend[Rest[x],#]&,List@@(x[[1]])],
		True,{x}]
	]/@triples
]


CollectRight[triples_List]:=
Module[{x,y},
	Table[Append[x,Total[Cases[triples,{Sequence@@x,y_}:>y]]],{x,DeleteDuplicates[Most/@triples]}]
]


ExpandRight[triples_List]:=
Module[{x},
	Function[x,Sequence@@Which[
		x[[-1]]===0,{},
		x[[-1,0]]===Plus,Map[Append[Most[x],#]&,List@@(x[[-1]])],
		True,{x}]
	]/@triples
]


(* ::Subsection::Closed:: *)
(*Adjungate operator definition*)


(* ::Text:: *)
(*Properties of the adjoint operator. Compatible with Prod[] and Mathematica's non-commutative multiplication.*)


adj[Prod[a_]]:=Prod[adj[a]]
adj[Prod[a__,b__]] := Prod[adj[Prod[b]],adj[Prod[a]]]
adj[Prod[]]:= Prod[]
adj[Times[a__,Prod[b___]]]:= a adj[Prod[b]]
adj[a___,b_Plus,c___]:=(adj[a,#,c]&/@b)
adj[adj[a__]] := a


adj[NonCommutativeMultiply[a_,b_]] := NonCommutativeMultiply[adj[b],adj[a]]
adj[1]:= 1
adj[-a_]:= -adj[a];
adj[Times[a_?CoeffQ,NonCommutativeMultiply[b___]]]:= a adj[NonCommutativeMultiply[b]]
adj[Times[a_?CoeffQ,b_]]:=a adj[b]


(* ::Subsection::Closed:: *)
(*Quiver*)


(* ::Text:: *)
(*Data structure*)


Quiver = {RepeatedNull[{__,__,__}]};


(* ::Text:: *)
(*Gives all Sources and Targets, respectively, of a certain lable of a quiver Q (not necessarily with unique lables).*)


Target[Alternatives[l_Symbol,l:adj[_]],Q:Quiver]:=
Module[{q},
	q = Cases[Q,{l,__,__}];
	If[Length[q]===0,{},q[[All,3]]]
]


Source[Alternatives[l_Symbol,l:adj[_]],Q:Quiver]:=
Module[{q},
	q = Cases[Q,{l,__,__}];
	If[Length[q]===0,{},q[[All,2]]]
]


(* ::Text:: *)
(*Returns the signature of a polynomial w.r.t. to a quiver Q (not necessarily with unique lables)*)


QSignature[l_List,Q:Quiver]:= Map[QSignature[#,Q]&,l]


QSignature[p_,Q:Quiver]:=
Module[{monomials,m,results,i,begin,end,comb,vertices},
	monomials = (MonomialList[ToProd[p]]/.c__*Prod[x___]->Prod[x])/.Prod->List;
	Catch[
	(*base case: only one mononmial*)
	If[Length[monomials] === 1,
			m = monomials[[1]];
			(*empty monomial = Prod[] \[Rule] is compatible; Return all empty paths*)
			If[Length[m]===0, 
				vertices = Sort[DeleteDuplicates[Flatten[Q[[All,2;;3]]]]];
				Throw[{Map[{#,#}&,vertices]}]];
			(*get all sources and targets of first lable*)
			begin = Cases[Q,{m[[-1]],_,_}][[All,2;;3]];
			For[i = Length[m]-1, i > 0, i--,
				(*get all sources and targets of lable i*)
				end = Cases[Q,{m[[i]],_,_}][[All,2;;3]];
				(*check if the two lables are compatible*)
				comb = {};
				Do[Map[If[i[[2]]===#[[1]],AppendTo[comb,{i[[1]],#[[2]]}]]&,end],{i,begin}];
				If[comb === {}, Throw[{}]];
				begin = comb;
			];
			Throw[begin],	
			(*usual case: split polynomial into monomials*)
			results = Map[QSignature[Prod@@#,Q]&,monomials];
			If[AllTrue[results,# === results[[1]]&],
				Throw[results[[1]]],
				Throw[{}]
			];
	]
	]
]


PlotQuiver[Q:Quiver]:=
	GraphPlot[Map[{#[[2]]->#[[3]],#[[1]]}&,Q],DirectedEdges->True,SelfLoopStyle->.2]


(* ::Subsection::Closed:: *)
(*Certify*)


Certify[assumptionsInput_List,claims_,Q:Quiver,OptionsPattern[{MaxIter->10,MaxDeg->Infinity,MultiLex->False,Info->False,Parallel->True,Sorted->False,Criterion->True}]]:=
 Module[{info,maxiter,reduced,vars,cofactors,G,assump,claim,certificate,rules,lc,toIgnore,toIgnoreOld,zeros,i,knowns,unknowns,t,assumptions,redCofactors,k,l,count,assumptionsRed},
	info = OptionValue[Info];
	maxiter = OptionValue[MaxIter];
	
	assumptions = ToProd/@assumptionsInput;
	
	(*check compatibility of the assumptions and the claims*)
	Do[
		If[QSignature[assump]==={},
			Print["The assumption ", assump , " is not compatible with the quiver."];
			Return[$Failed]
		]
	,{assump,assumptions}
	];
	
	If[Head[claims] === List,
		Do[
			If[QSignature[claim]==={},
				Print["The claim ", claim , " is not compatible with the quiver."];
				Return[$Failed]
			]
		,{claim,claims}
		],
		If[QSignature[claims,Q] === {},
		Print["The claim is not compatible with the quiver."]; Return[$Failed]]
	];
	
	(*set up the ring*)
	Print["Using the following monomial ordering:"];
	If[OptionValue[MultiLex],
		knowns = DeleteDuplicates[If[Head[claims]===List,
					Cases[Q[[All,1]],Alternatives@@Flatten[Map[#/.Alternatives[Plus,Times,NonCommutativeMultiply]->List&,claims]]],
					Cases[Q[[All,1]],Alternatives@@Flatten[claims/.Alternatives[Plus,Times,NonCommutativeMultiply]->List]]
					]];
		unknowns = DeleteDuplicates[Cases[Q[[All,1]],var_/;!MemberQ[knowns,var]]];
		SetUpRing[knowns,unknowns],
		SetUpRing[DeleteDuplicates[Q[[All,1]]]]
	];

	(*make ideal monic*)
	lc = MakeMonic[assumptions];
	
	(*interreduce the generators*)
	{assumptionsRed,redCofactors} = Interreduce[assumptions,InputProd->True];
	If[info, Print["\nInterreduced the input from ", Length[assumptionsInput], " polynomials to ", Length[assumptionsRed], ".\n"]];
	
	(*compute the Groebner basis and reduce the claims*)
	If[info, Print["Computing a (partial) Groebner basis and reducing the claim...\n"]];
	(*do computation iteratively*)
	cofactors = {};
	If[Head[claims]===List,zeros = ConstantArray[0,Length[claims]],zeros=0];
	toIgnore = Length[assumptionsRed];
	i = 1;
	If[info,Print["Starting iteration ", i++ ,"...\n"]];
	G = Groebner[cofactors,assumptionsRed,1,MaxDeg->OptionValue[MaxDeg],Info->OptionValue[Info],Parallel->OptionValue[Parallel],Sorted->OptionValue[Sorted],Criterion->OptionValue[Criterion],OutputProd->True,Rewrite->False];
	reduced = ReducedForm[vars,G,claims];
	While[reduced =!= zeros && i <= maxiter,
		toIgnoreOld = Length[G];
		If[info,Print["Starting iteration ", i++ ,"...\n"]];
		G = Groebner[cofactors,G,1,Ignore->toIgnore,MaxDeg->OptionValue[MaxDeg],Info->OptionValue[Info],Parallel->OptionValue[Parallel],Sorted->OptionValue[Sorted],Criterion->OptionValue[Criterion],OutputProd->True,Rewrite->False];
		toIgnore = toIgnoreOld;
		count = 0;
		While[reduced =!= zeros && count < 4,
			reduced = ReducedForm[vars,RandomSample[G],claims];
			count++;
		];
	];
	If[info, Print["Rewriting the cofactors has started..."]];
	t = AbsoluteTiming[RewriteGroebner[cofactors,Info->OptionValue[Info],OutputProd->True]][[1]];
	If[info, Print["Rewriting the cofactors took in total ", t]];
	
	(*rewrite the linear combination*)
	If[info,
		Print["\nRewriting the linear combination in terms of the assumptions has started...\n"]];
	certificate = Rewrite[vars,cofactors,InputProd->True];
	
	(*rewrite in terms of assumptions and not of the interreduced assumptions*)
	rules = Table[{k_,assumptionsRed[[i]],l_}->
		Sequence@@Table[{Prod[k,j[[1]]],j[[2]],Prod[j[[3]],l]},{j,redCofactors[[i]]}],{i,Length[redCofactors]}];
	certificate = certificate/.rules;
	(*take care of leading coefficients in the certificate*)
	rules = MapIndexed[{a_,#1,b_}->{a/lc[[First[#2]]],Expand[lc[[First[#2]]]*#1],b}&,assumptions];
	certificate = certificate/.rules;
	(*convert back to NonCommutativeMultiply*)
	certificate = If[Head[claims]===List,
		ToNonCommutativeMultiply[Map[ExpandLeft[CollectLeft[#]]&,certificate]],
		ToNonCommutativeMultiply[ExpandLeft[CollectLeft[certificate]]]
	];
	(*return the reduced claims and the linear combinations*)
	If[info,
		If[reduced === zeros,
			Print["Done! All claims were successfully reduced to 0."],
			Print["Done! Not all claims could be reduced to 0."]
		]
	];
	{sigAssump,sigClaim,reduced,certificate}
]


(* ::Subsection::Closed:: *)
(*End*)


Copyright[a_String,b___String]:=Print[StringJoin[Prepend[{"\n",#}&/@{b},a]]]


Copyright[
    "Package OperatorGB version 1.0.1",
    "Copyright 2019, Institute of Algebra, JKU",
    "written by Clemens Hofstadler"];


(*End[]*)


EndPackage[]


?OperatorGB
