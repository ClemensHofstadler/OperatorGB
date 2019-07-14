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


ToProd[poly_]:= Prod[(poly//.NonCommutativeMultiply->Prod)]


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
(*Data Type*)


ReductionSystem={RepeatedNull[List[List[___],Function[___]]]};


(* ::Subsection:: *)
(*Ambiguities*)


(* ::Text:: *)
(*Reducible Words*)


ExtractReducibleWords[sys:ReductionSystem]:=
	MapIndexed[{#1[[1]],Sequence@@#2}&,sys]


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
Module[{k},
	Reap[For[k=1,k<Length[v]&&k<Length[w],k++,
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


DeleteRedundant[amb_List]:=
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




(* ::Subsection::Closed:: *)
(*S-Polynomials*)


(* ::Text:: *)
(*SPoly[amb,fi,fj] computes the S-polynomial corresponding to the ambiguity amb, which comes from the two reduction rules fi and fj. Additionally, the linear combination how the S-polynomial was computed from fi and fj is returned in a list. *)


SPoly[amb:(Alternatives[Overlap,Inclusion][_List,_List,_List,_List]),fi:{List[___],Function[___]},fj:{List[___],Function[___]}]:=
Module[{A,C},

	C = Prod@@amb[[2]];
	A = Prod@@amb[[3]];
	If[amb[[0]]=== Overlap,
			(*Overlap[ABC,C,A]*)
			{Prod[fi[[2]][Sequence@@fi[[1]]],C] - Prod[A,fj[[2]][Sequence@@fj[[1]]]],
				{{A,ToPoly[fj],Prod[]},{-Prod[],ToPoly[fi],C}}},
			(*Inclusion[CBA,C,A]*)
			{fi[[2]][Sequence@@fi[[1]]] - Prod[C,fj[[2]][Sequence@@fj[[1]]],A],
			{{C,ToPoly[fj],A},{-Prod[],ToPoly[fi],Prod[]}}}
	]
]


(* ::Text:: *)
(*Same as Poly but without returning the linear combination.*)


SPoly2[amb:(Alternatives[Overlap,Inclusion][_List,_List,_List,_List]),fi:{List[___],Function[___]},fj:{List[___],Function[___]}]:=
Module[{A,C},

	C = Sequence@@amb[[2]];
	A = Sequence@@amb[[3]];

	If[amb[[0]]=== Overlap,
			(*Overlap[ABC,C,A]*)
			Prod[fi[[2]][Sequence@@fi[[1]]],C] - Prod[A,fj[[2]][Sequence@@fj[[1]]]],
			(*Inclusion[CBA,C,A]*)
			Prod[fi[[2]][Sequence@@fi[[1]]]] - Prod[C,fj[[2]][Sequence@@fj[[1]]],A]
	]
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
Module[{x,y,count,spol,lt,info,p,h,G,r,t1,t2,lists,rules,sorted,oldlength,parallel,hrule,syslength,pos,incl,possible,maxdeg,outputProd,rulesCrit,linearComb},
info = OptionValue[Info];
sorted = OptionValue[Sorted];
parallel = OptionValue[Parallel];
maxdeg = OptionValue[MaxDeg];
outputProd = OptionValue[OutputProd];

If[Head[cofactors]=!=List,cofactors={}];

G = CreateRedSys[ideal];
oldlength = syslength = Length[G];
t1 = 0; t2 = 0; count = 0;
If[info,Print["G has ", Length[G]," elements in the beginning."],Print[]];

spol = DeleteDuplicates[CheckResolvability[G,OptionValue[Ignore],Criterion->OptionValue[Criterion],MaxDeg->maxdeg,Info->info,Sorted->sorted,Parallel->parallel],SameQ[#1[[1]],#2[[1]]]&];
rules = ExtractRules[G];

While[Length[spol] > 0 && count < maxiter,
	t1 = AbsoluteTiming[
	Monitor[While[Length[spol]\[NonBreakingSpace]> 0,
		
		(*pick next S-polynomial*)
		p = First[spol];
		spol = Drop[spol,1];
		
		(*reduce it*)
		r = Reap[p[[1]]//.rules]; 
		h = r[[1]];
		If[Length[r[[2]]] > 0,
			p[[2]]\[NonBreakingSpace]= Join[p[[2]],r[[2,1]]]
		];
		lt = LeadingTerm[h];
		If[lt[[1]] =!= 1, 
			h = Expand[1/lt[[1]]*h]; p[[2]] = (ReplacePart[#,1 -> 1/lt[[1]]*#[[1]]]&/@ p[[2]])];
		If[h =!= 0,
			hrule = CreateRedSys[h];
			AppendTo[G,hrule];
			AppendTo[cofactors,{h,p[[2]]}]; 
			AppendTo[rules,Sequence@@ExtractRules[{hrule}]];
		];
	];,Length[spol]];][[1]];

	If[info, Print["The reduction took ", t1]];
	count = count + 1;
	If[info,Print["Iteration ",count, " finished. G has now ", Length[G]," elements\n"]];
	If[count < maxiter, 
		spol = DeleteDuplicates[CheckResolvability[G,oldlength,Criterion->OptionValue[Criterion],MaxDeg->maxdeg,Info->info,Sorted->sorted,Parallel->parallel],SameQ[#1[[1]],#2[[1]]]&];
		oldlength = Length[G];
	];
];
If[OptionValue[Rewrite],
	If[info, Print["Rewriting the cofactors has started."]];
	t2 = AbsoluteTiming[RewriteGroebner[cofactors,Info->info,OutputProd->outputProd];][[1]];
	If[info, Print["Rewriting the cofactors took in total ", t2]];
];
If[outputProd,
	ToPoly[G],
	Map[ToNonCommutativeMultiply,ToPoly[G]]
]
]


(* ::Text:: *)
(*CheckResolvability[sys,Info->False,Parallel->True,Sorted->False] returns all S-polynomials from the reduction system sys which can not be reduced to zero. Additionally, *)
(*for each S-polynomial a list containing the linear combination how the S-polynomial was generated from the elements of sys is returned.*)
(*For a description of the OptionPatterns see the documentation of the Groebner method.*)


CheckResolvability[sys:ReductionSystem,oldlength:_?IntegerQ:0,OptionsPattern[{Criterion->False,MaxDeg->Infinity,Info->False,Parallel->True,Sorted->False}]]:=
Module[{amb,spol,info,t1,t2,rules,sorted,maxdeg,parallel,words,x,y,r},

	info = OptionValue[Info];
	sorted = OptionValue[Sorted];
	parallel = OptionValue[Parallel];
	maxdeg = OptionValue[MaxDeg];

	(*generate ambiguities*)
	words = ExtractReducibleWords[sys];
	t1 = AbsoluteTiming[
		amb = GenerateAmbiguities[words[[;;oldlength]],words[[oldlength+1;;]],maxdeg,Parallel->parallel];
	][[1]];
	If[info,Print[Length[amb]," ambiguities in total (computation took ",t1, ")"]];
	
	If[OptionValue[Criterion],
		amb = DeleteRedundant[amb]
	];
	
	If[sorted,amb = Sort[amb]];
	
	
	(*generate S-polynomials*)
	t2 = AbsoluteTiming[
	x = AbsoluteTiming[
	If[parallel,
		spol = DeleteCases[ParallelMap[SPoly[#,sys[[#[[4,1]]]],sys[[#[[4,2]]]]]&,amb,DistributedContexts->Automatic,Method->"CoarsestGrained"],{0,___}],
		spol = DeleteCases[Map[SPoly[#,sys[[#[[4,1]]]],sys[[#[[4,2]]]]]&,amb],{0,___}]
	];		
	][[1]];
	If[info,Print["Generating S-polys: ",x]];
	(*reduce S-polynomials*)
	rules = ExtractRules[sys];
	
	(*parallelizing this makes it only slower*)
	y = AbsoluteTiming[
	spol = DeleteCases[Map[(r = Reap[#[[1]]//.rules]; If[r[[1]]=!= 0,
											If[Length[r[[2]]] > 0,
												{r[[1]],Join[#[[2]],r[[2,1]]]},
												{r[[1]],#[[2]]}
												],
											{}
											])&,spol],{}];
	][[1]];
	If[info,Print["Reducing S-polys: ",y]];
	][[1]];
	If[info, Print[Length[spol]," different S-polynomials did not reduce to 0 (computation took ",t2,")"]];
	spol
]


(* ::Text:: *)
(*ExtractRules[vars,sys] generates a set of reduction rules out of the reduction system sys with the special feature that the cofactors of the reduction will be sowed whenever such a rule is applied. They can then be reaped using the Reap command.*)
(*Sys has to consist of pairs {word,func}, which can be obtained using the method CreateRedSys.*)


ExtractRules[sys:ReductionSystem]:=
Module[{a,b,c,coeff,h,i,j,p,terms,q},
	a=Unique[];b=Unique[];
	coeff = Unique[];
	terms = sys[[All,1]];
	Table[
		p=Prod[Pattern[Evaluate[a],BlankNullSequence[]],terms[[i]],Pattern[Evaluate[b],BlankNullSequence[]]];
		q = Expand[Evaluate[coeff]*Prod[Evaluate[a],sys[[i,2]][Sequence@@terms[[i]]],Evaluate[b]]];
		With[{x ={-Evaluate[coeff]*Prod[Evaluate[a]],ToPoly[sys[[i]]],Prod[Evaluate[b]]},y=q},Alternatives[Pattern[Evaluate[coeff],BlankNullSequence[]]*p,p]:> Hold[Sow[x];y]]
	,{i,Length[sys]}]//ReleaseHold
];


SetAttributes[RewriteGroebner,HoldFirst]

RewriteGroebner[cofactors_,OptionsPattern[{Info->False,OutputProd->False}]]:=
Module[{a,b,i,j,count,rules,info,occurring},
	info = OptionValue[Info];
	rules = {};
	count = Length[cofactors];
	
	Monitor[Do[
		(*only choose those rules which will actually be used*)
		occurring= Position[cofactors[[;;i,1]],Alternatives@@DeleteDuplicates[cofactors[[i+1,2,All,2]]],{1}];

		rules = Map[{a__,#[[1]],b__} -> Sequence@@Table[{Prod[Evaluate[a],#[[2,j,1]]],#[[2,j,2]],Prod[#[[2,j,3]],Evaluate[b]]},{j,Length[#[[2]]]}]&,Extract[cofactors,occurring]];
		
		(*reduce the cofactorlist*)
		cofactors[[i+1,2]] = cofactors[[i+1,2]]/.rules//CollectLeft//ExpandLeft;

		count = count - 1;
	,{i,Length[cofactors]-1}];
	,count];
	If[!OptionValue[OutputProd],
		cofactors = Map[Map[ToNonCommutativeMultiply,#]&,cofactors];
	];
]


(* ::Subsection::Closed:: *)
(*Find Cofactors*)


(* ::Text:: *)
(*ReducedForm[cofactors,G,exp] can be used to reduce the expression exp with the elements of G. The linear combination of these reduction steps is saved in the list cofactors. *)
(*Exp can also be a list of expressions. *)


SetAttributes[ReducedForm,HoldFirst]

ReducedForm[cofactors_,G_,exp_]:=
Module[{t,lists,rules,sys},
	cofactors = {};
	sys = CreateRedSys[G];
	rules = ExtractRules[sys];
	t = Reap[ToProd[exp]//.rules];
	If[Length[t[[2]]]>0,
		cofactors = ReplacePart[#,1->-#[[1]]]&/@t[[2,1]];
		cofactors = Map[ToNonCommutativeMultiply,cofactors]
	];
	ToNonCommutativeMultiply[t[[1]]]
]


ReducedForm[cofactors_,G_,exp:_?ListQ]:=
Module[{i,j,k,t,lists,rules,sys},
	cofactors = Table[{},{i,Length[exp]}];
	sys = CreateRedSys[G];
	rules = ExtractRules[sys];
	Table[
		t = Reap[ToProd[exp[[j]]]//.rules];
		If[Length[t[[2]]]>0,
			cofactors[[j]] = ReplacePart[#,1->-#[[1]]]&/@t[[2,1]];
			cofactors[[j]] = Map[ToNonCommutativeMultiply,cofactors[[j]]];
		];
		ToNonCommutativeMultiply[t[[1]]],{j,Length[exp]}
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
Module[{count,spol,p,h,G,lt,info,t,rules,criterion,oldlength,maxdeg,multiples,incl,hrule,pos,sorted,parallel,syslength},

info = OptionValue[Info];
criterion = OptionValue[Criterion];
sorted = OptionValue[Sorted];
parallel = OptionValue[Parallel];
maxdeg = OptionValue[MaxDeg];

G = CreateRedSys[ideal];
syslength = oldlength = Length[G];
If[info,Print["G has ", Length[G]," elements in the beginning."];Print[]];
count = 0; t = 0;

spol = DeleteDuplicates[CheckResolvability2[G,OptionValue[Ignore],MaxDeg->maxdeg,Info->info,Sorted->sorted,Parallel->parallel]];
rules = ExtractRules2[G];

While[Length[spol] > 0 && count < maxiter,
	t = AbsoluteTiming[Monitor[While[Length[spol] > 0,
		p = First[spol];
		spol = Drop[spol,1];
		h = p//.rules;
		lt = LeadingTerm[h];
		If[lt[[1]] =!= 1, h = Expand[1/lt[[1]]*h]];
		If[h =!= 0, 
			hrule = CreateRedSys[h];
			If[criterion,
				pos = Position[G[[syslength+1;;]],{{___,Sequence@@lt[[2]],___},__}];
				incl = Flatten[DeleteCases[Map[Inclusion[{G[[#,1]],#},{lt[[2]],Length[G]+1}]&,Flatten[pos]],{}]];
				pos = Cases[Map[{SPoly2[#,G[[#[[4,1]]]],hrule],#[[4,1]]}&,incl],{0,___}];
				pos = Partition[pos[[All,2]],1];
				G = Delete[G,pos]; rules = Delete[rules,pos];
			];
			AppendTo[G,hrule]; AppendTo[rules,Sequence@@ExtractRules2[{G[[-1]]}]]];
	];,Length[spol]];][[1]];
	If[info, Print["The reduction took ", t]];
	count = count + 1;
	If[info,Print["Iteration ",count, " finished. G has now ", Length[G]," elements"];Print[]];
	If[count < maxiter,
		spol = DeleteDuplicates[CheckResolvability2[G,oldlength,MaxDeg->maxdeg,Info->info,Sorted->sorted,Parallel->parallel]];
		oldlength = Length[G];
	];];
Map[ToNonCommutativeMultiply,ToPoly[G]]
]


CheckResolvability2[sys:ReductionSystem,oldlength:_?IntegerQ:0,OptionsPattern[{MaxDeg->Infinity,Info->False,Sorted->False,Parallel->True}]]:=
Module[{amb,spol,info,t1,t2,t3,lists,rules,words,sorted,parallel,maxdeg},
	info = OptionValue[Info];
	sorted = OptionValue[Sorted];
	parallel = OptionValue[Parallel];
	maxdeg = OptionValue[MaxDeg];

	words = ExtractReducibleWords[sys];
	rules = ExtractRules2[sys];

	(*generate ambiguities*)
	t1 = AbsoluteTiming[
		amb = GenerateAmbiguities[words[[;;oldlength]],words[[oldlength+1;;]],maxdeg,Parallel->True]
	][[1]];
	If[info,Print[Length[amb]," ambiguities in total (computation took ", t1, ")"]];
	
	If[sorted,amb = Sort[amb]];
	
	(*generate and reduce S-polynomials*)
	t2 = AbsoluteTiming[
	If[parallel && Length[amb] > 300,
		spol = DeleteCases[ParallelMap[SPoly2[#,sys[[#[[4,1]]]],sys[[#[[4,2]]]]]//.rules &,amb,DistributedContexts->Automatic,Method->"ItemsPerEvaluation" -> 1000],0],
		spol = DeleteCases[Map[SPoly2[#,sys[[#[[4,1]]]],sys[[#[[4,2]]]]]//.rules &,amb],0]];
	][[1]];

	If[info, Print[Length[spol]," different S-polynomials did not reduce to 0 (computation took ",t2,")"]];
	spol
]


(* ::Text:: *)
(*Parallelising this function does not make sense. It is faster when executed sequential.*)


ExtractRules2[sys:ReductionSystem]:=
Module[{a,b,h,i,j,m,p,terms},
	a=Unique[];b=Unique[];
	terms = sys[[All,1]];
	Table[Prod[Pattern[Evaluate[a],BlankNullSequence[]],terms[[i]],Pattern[Evaluate[b],BlankNullSequence[]]]->
		Prod[a,sys[[i,2]][Sequence@@terms[[i]]],b],{i,Length[terms]}]
]


ApplyRules[expr_,G_]:= Module[
{sys},
	sys = CreateRedSys[G];
	ToNonCommutativeMultiply[ToProd[expr]//.ExtractRules[sys]]
]


(* ::Subsection::Closed:: *)
(*Additional stuff*)


NormalizePoly[poly_]:= 
	Expand[1/LeadingTerm[poly][[1]]*poly];


Remainder[poly_, lt:List[__,List[___]]]:=
	poly - lt[[1]]*Prod@@lt[[2]]


CreateRedSys[polies_List]:=
	Map[CreateRedSys,polies]


CreateRedSys[p_]:=
Module[{m,lt,poly,i},
	poly = ToProd[p];
	lt = LeadingTerm[poly];
	m=Table[Unique[],{Length[lt[[2]]]}];
	{lt[[2]], Function@@{Evaluate/@m , Expand[-1/lt[[1]]*Remainder[poly,lt]]}//.Table[lt[[2]][[i]] -> m[[i]],{i,Length[m]}] }
]


ToPoly[poly:{List[__],Function[___]}]:= Prod@@poly[[1]] - poly[[2]][Sequence@@poly[[1]]]


ToPoly[sys:ReductionSystem]:= ToPoly/@sys;


Rewrite[spolfactors:List[RepeatedNull[List[RepeatedNull[{__,__,__}]]]],cofactors_,OptionsPattern[InputProd->False]]:= 
	Map[Rewrite[#,cofactors,InputProd->OptionValue[InputProd]]&,spolfactors]


Rewrite[spolfactor_, cofactor_,OptionsPattern[InputProd->False]]:=
Module[{a,b,i,j,rules,occurring,result,spolfactors,cofactors},
	spolfactors = Map[{ToProd[#[[1]]],ToProd[#[[2]]],ToProd[#[[3]]]}&,spolfactor];
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
	Map[ToNonCommutativeMultiply, result//CollectLeft//ExpandLeft]
]


MultiplyOut[cofactors_List]:=Expand[ToNonCommutativeMultiply[Total[Map[ToProd,cofactors]]]]


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
adj[Times[a__,NonCommutativeMultiply[b___]]]:= a adj[NonCommutativeMultiply[b]]


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
Module[{monomials,m,results,s,t,i,begin,end,comb,poly,vertices},
	poly = ToProd[p];
	monomials = (MonomialList[poly]/.c__*Prod[x___]->Prod[x])/.Prod->List;
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


Certify[assumptions_List,claims_,Q:Quiver,OptionsPattern[{MaxIter->10,MaxDeg->Infinity,MultiLex->False,Info->False,Parallel->True,Sorted->False}]]:=
 Module[{info,maxiter,reduced,vars,cofactors,G,sigAssump,sigClaim,certificate,rules,lc,toIgnore,toIgnoreOld,zeros,i,knowns,unknowns,t},
	info = OptionValue[Info];
	maxiter = OptionValue[MaxIter];
	
	(*check compatibility of the assumptions and the claims*)
	sigAssump = Map[QSignature[#,Q]&,assumptions];
	If[MemberQ[sigAssump,{}],
		Print["One of the assumptions is not compatible with the quiver."]; Return[$Failed]];
	If[Head[claims] === List,
		sigClaim = Map[QSignature[#,Q]&,claims];
		If[MemberQ[sigClaim,{}], Print["One of the claims is not compatible with the quiver."]; Return[$Failed]],
		sigClaim = QSignature[claims,Q];
		If[sigClaim === {},
		Print["The claim is not compatible with the quiver."]; Return[$Failed]]
	];
	
	(*set up the ring*)
	If[info,
		Print["Using the following monomial ordering:"]];
	If[OptionValue[MultiLex],
		knowns = DeleteDuplicates[If[Head[claims]===List,
					Cases[Q[[All,1]],Alternatives@@Flatten[Map[#/.Alternatives[Plus,Times,NonCommutativeMultiply]->List&,claims]]],
					Cases[Q[[All,1]],Alternatives@@Flatten[claims/.Alternatives[Plus,Times,NonCommutativeMultiply]->List]]
					]];
		unknowns = DeleteDuplicates[Cases[Q[[All,1]],var_/;!MemberQ[knowns,var]]];
		SetUpRing[knowns,unknowns],
		SetUpRing[DeleteDuplicates[Q[[All,1]]]]
	];
	
	(*compute the Groebner basis and reduce the claims*)
	If[info, Print["\n","Computing a (partial) Groebner basis and reducing the claim...\n"]];
	(*do computation iteratively*)
	cofactors = {};
	If[Head[claims]===List,zeros = ConstantArray[0,Length[claims]],zeros=0];
	toIgnore = Length[ideal];
	i = 1;
	If[info,Print["Starting iteration ", i++ ,"...\n"]];
	G = Groebner[cofactors,assumptions,1,MaxDeg->OptionValue[MaxDeg],Info->OptionValue[Info],Parallel->OptionValue[Parallel],Sorted->OptionValue[Sorted],OutputProd->True,Rewrite->False];
	reduced = ReducedForm[vars,G,claims];
	While[reduced =!= zeros && i <= maxiter,
		toIgnoreOld = Length[G];
		If[info,Print["Starting iteration ", i++ ,"...\n"]];
		G = Groebner[cofactors,G,1,Ignore->toIgnore,MaxDeg->OptionValue[MaxDeg],Info->OptionValue[Info],Parallel->OptionValue[Parallel],Sorted->OptionValue[Sorted],OutputProd->True,Rewrite->False];
		toIgnore = toIgnoreOld;
		reduced = ReducedForm[vars,G,claims];
	];
	If[info, Print["Rewriting the cofactors has started..."]];
	t = AbsoluteTiming[RewriteGroebner[cofactors,Info->OptionValue[Info],OutputProd->True]][[1]];
	If[info, Print["Rewriting the cofactors took in total ", t]];
	
	(*rewrite the linear combination*)
	If[OptionValue[Info],
		Print["\nRewriting the linear combination in terms of the assumptions has started..."]];
	certificate = Rewrite[vars,cofactors,InputProd->True];
	(*take care of leading coefficients in the certificate*)
	rules = Map[(lc = LeadingTerm[#][[1]];{a_,#/lc,b_}->{a/lc,#,b})&,assumptions];
	If[Head[claims]===List,
		certificate = Map[#/.rules&, certificate],
		certificate = certificate/.rules
	];
	
	(*return the reduced claims and the linear combinations*)
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
