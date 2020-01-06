(* ::Package:: *)

BeginPackage["OperatorGB`"]


Global`OperatorGB::usage="";


Clear[
	CoeffQ,Prod,ToProd,ToNonCommutativeMultiply,
	SetUpRing,
	WordOrder,varSets,
	LeadingTerm,DegLex,WeightedDegLex,MultiLex,Weight,SortedQ,
	ExtractReducibleWords,GenerateAmbiguities,
	Groebner,
	F4,
	ReducedForm,
	GroebnerWithoutCofactors,ApplyRules,
	CreateRedSys,ToPoly,Rewrite,MultiplyOut,Interreduce,
	CollectLeft,CollectRight,ExpandLeft,ExpandRight,
	adj,Pinv,AddAdj,
	Quiver,QSignature,PlotQuiver,CompatibleQ,UniformlyCompatibleQ,QOrderCompatibleQ,QConsequenceQ,
	Certify,MaxIter,MaxDeg,Info,Parallel,Sorted,Criterion
]


(*Non-commutative multiplication*)
CoeffQ::usage="CoeffQ[k_] should check whether k is in the base ring of the polynomial ring.";
Prod::usage="Prod[m1,m2,...] represents the non-commutative multiplication.";
ToProd::usage="Convert a polynomial from the built in non-commutative multiplication to the Prod data structure";
ToNonCommutativeMultiply::usage="Convert a polynomial from the Prod data structure to the built in non-commutative multiplication";


(*SetUpRing*)
SetUpRing::usage="SetUpRing defines the non-commutative polynomial ring and a monomial order.

This method can be called with one list X = {x1,x2,...,xn} containing all indeterminates the user wants to work with. This then sets up the non-commutative polynomial ring \!\(\*TemplateBox[{},\n\"Rationals\"]\)<X> of non-commutative polynomials in the variables x1,x2,...,xn with rational numbers as coefficients. By default, this also sets up a degree lexicographic order where the indeterminates are ordered in increasing order according to their appearance in the list X, i.e. x1 < x2 < ... < xn.

The user has the option to change the monomial order. The second pre-defined order is a weighted degree lexicographic order. To switch to this order, the user has execute the command SortedQ := WeightedDegLex. Then, additionally, a list of weights {w1,w2,...} named Weight has to be provided. Type ?Weight for more information.

To define an individual order, the user can provide a binary function SortedQ[a,b] working on the set of all words that can be built from the alphabet of indeterminates in X. Type ?SortedQ for more information.

It is also possible to set up a multigraded lexicographic order. To do this, the user has to call SetUpRing with two lists X = {x1,x2,...,xn} and Y = {y1,y2,...,ym} as input. This defines the non-commutative polynomial ring \!\(\*TemplateBox[{},\n\"Rationals\"]\)<X,Y> of non-commutative polynomials in the variables x1,x2,...,xn,y1,y2,...,ym with rational numbers as coefficients. In this case, by default a multigraded lexicographic order is defined, where two monomials M1 and M2 are first compared by their degree only in indeterminates from Y, then by their degree only in indeterminates from X and finally, to break ties, a degree lexicographic order x1 < ... < xn < y1 < ... < ym is used. The multigraded lexicographic order is denoted by x1 < x2 < ... < xn << y1 < y2 < ... < ym.
"


(*Ordering*)
LeadingTerm::usage="Leading term of the polynomial w.r.t. the specified monomial ordering."
DegLex::usage="Degree Lexicographic order"
WeightedDegLex::usage="Weighted Degree Lexicographic order"
MultiLex::usage="Multigraded Lexicographic order"
Weight::usage="A list {w1,w2,...} defining weights on the variables for the weighted degree lexicographic order. The weight w1 corresponds to the first variable in the input list of SetUpRing, w2 to the second one and so on."
SortedQ::usage="SortedQ[M1,M2] is a binary function working on the set of all words that can be built from the alphabet of variables given in SetUpRing. It decides, when given two words M1 and M2 as input, which of them is larger. SortedQ[M1,M2] returns True if M1 \[LessEqual] M2 and False otherwise.

M1 and M2 have to be given in form of lists containig only elements from the list(s) X (and Y) that where used in the call of the method SetUpRing. For example, if X = {x,y,z} then M1 could be of the form {x,x,y,z} and M2 could be {y,z,x,y}.
"


(*ambiguities*)
ExtractReducibleWords::usage="Preprocessing before GenerateAmbiguities."
GenerateAmbiguities::usage="GenerateAmbiguities[words] computes all ambigutites among all words in the set 'words'."


(*Groebner basis*)
Groebner::usage="Groebner[cofactors_,ideal_, maxiter:_?IntegerQ:10, OptionsPattern[{Criterion->True,Ignore->0,MaxDeg->Infinity,Info\[Rule]True,Parallel->True,Sorted->True}]] executes at most maxiter iterations of the Buchberger algorithm to compute
a (partial) Groebner basis of an ideal. Additionally, for every new element in the Groebner basis a list of cofactors is saved in the list cofactors forming a linear combination of the new element. For further information concerning the OptionPatterns
please see the documentation or the source code."


(*F4)*)
F4::usage="F4[cofactors_,ideal_, maxiter:_?IntegerQ:10, OptionsPattern[{N->100,Criterion->True,Ignore->0,MaxDeg->Infinity,Info\[Rule]True,Parallel->True,Sorted->True}]] executes at most maxiter iterations of Faugere's F4 algorithm to compute
a (partial) Groebner basis of an ideal. Additionally, for every new element in the Groebner basis a list of cofactors is saved in the list cofactors forming a linear combination of the new element. For further information concerning the OptionPatterns
please see the source code."


(*Find cofactors*)
ReducedForm::usage="ReducedForm[cofactors,G,exp] reduces the expression exp by the elements of G and saves the cofactors of the reduction process in the list cofactors. The argument exp can also
be a list of expressions, then all expressions are reduced."


(*Groebner basis without cofactors*)
GroebnerWithoutCofactors::usage="GroebnerWithoutCofactors[ideal_,maxiter:_?IntegerQ:10,OptionsPattern[{MaxDeg->Infinity,Info->False,Parallel->True,Sorted->True,Criterion->True}]] executes at most maxiter iterations 
of the Buchberger algorithm to compute a (partial) Groebner basis of an ideal. For further information concerning the OptionPatterns please see the documentation or the source code."
ApplyRules::usage="ApplyRules[exp,G] reduces the expression exp using polynomials from the set G."


(*Additional stuff*)
CreateRedSys::usage="Converts polynomials into the data structure needed for the Groebner algorithm."
ToPoly::usage="Converts an element of a reduction system back into a polynomial."
Rewrite::usage="Rewrite[vars, cofactors] rewrites a linear combination, which is safed in vars, with the elements from cofactors."
MultiplyOut::usage="To multiply out a list of cofactors given in terms of the built in non-commutative multiplication."
Interreduce::usage="Interreduce[ideal_] interreduces the polynomials in 'ideal'."


(*Collect cofactors*)
CollectLeft::usage="Tries to collect triples in a cofactor representation having the same left cofactors."
CollectRight::usage="Tries to collect triples in a cofactor representation having the same right cofactors."
ExpandLeft::usage=""
ExpandRight::usage=""


(*Adjungate operator definition*)
adj::usage="adj[A] represents the adjoint of the operator A"
Pinv::usage="Gives the 4 Moore-Penrose equations."
AddAdj::usage="Adds the adjoint statements to a given list of statements."


(*Quiver*)
Quiver::usage="Data structure of a Quiver"
QSignature::usage="QSignature[poly,Q] returns the signature of the polynomial poly w.r.t. the quiver Q (not necessarily with unique labels)"
PlotQuiver::usage="Plot a quiver Q"
CompatibleQ::usage="Tests whether a polynomial is compatible with a quiver."
UniformlyCompatibleQ::usage="Tests whether a polynomial is uniformly compatible with a quiver."
TrivialQuiver::usage="Returns the trivial quiver containing all variables of the given input."
QOrderCompatibleQ::usage="Tests whether a polynomial is Q-order-compatible with a quiver and the order defined by SetUpRing."
QConsequenceQ::usage="Tests whether a certificate is a Q-consequence of a set of polynomials and a quiver."


(*Certify*)
Certify::usage="Certifies whether a certain claim is a consequence of some assumptions via Groebner 
basis computations. Additionally, compatibility with a given quiver is checked."


Begin["`Private`"]


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


SetUpRing[vars_List,OptionsPattern[{Info->True}]]:= (
WordOrder = vars;
SortedQ := DegLex;
If[OptionValue[Info],
	Print[Sequence@@Map[ToString[#,StandardForm]<>" < "&,vars[[;;-2]]] <> ToString[vars[[-1]],StandardForm]];
];
);


SetUpRing[S1_List,SRest:_List..,OptionsPattern[{Info->True}]]:= Module[{string},
	WordOrder = Flatten[Join[S1,SRest]]; 
	varSets = {S1,SRest};
	SortedQ := MultiLex;
	If[OptionValue[Info],
		string = "";
		Do[
			string = string <> Sequence@@Map[ToString[#,StandardForm]<>" < "&,S[[;;-2]]] <> ToString[S[[-1]],StandardForm];
			string = string <> " << ";
		,{S,varSets}];
		Print[StringTake[string,{1,-4}]];
	];
	varSets = Reverse[varSets];
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


MultiLex[a_List,b_List]:= Module[{Va,Vb},
Catch[
	Do[
		Va = Count[a,Alternatives@@S];
		Vb = Count[b,Alternatives@@S];
		If[Va < Vb, Throw[True]];
		If[Va > Vb, Throw[False]];
	,{S,varSets}];
	Throw[DegLex[a,b]];
]
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
(*Hence, this function determines the monomial order. By defining an own SortedQ function the user can implement his own monomial order. By default the Degree Lexicographic order is set. It is also possible to switch to a weighted Degree Lexicographic order by setting SortedQ := WeightedDegLex. Then the user also has to define a list Weight = {w1,w2,...} where w1 defines the weight for the first module in WordOrder, w2 defines the weight for the second module in WordOrder and so on.*)


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
			Sow[Overlap[Join[v,Drop[w,k]],Drop[v,-k],Drop[w,k],{i,j}]]
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
(*Chain Criterion*)


DeleteRedundant[amb_List,lt_List,OptionsPattern[{Info->False}]]:= 
Module[{result,t,i,j,possibilities,V},
	t = AbsoluteTiming[
	i = 0;
	result = amb;
	Do[
		V = Sequence@@lt[[i]];
		possibilities = With[{idx = i},Alternatives@@{
			_[{___,V,___},_,_,ij_] /; idx < Min[ij], 
			Overlap[_,{___,V,___},_,_],
			Overlap[_,_,{___,V,___},_],
			Overlap[{_,___,V,___,_},_,_,_],
		     Inclusion[_,{___,V,___},_,{i_,j_}]/; i < idx,
		    Inclusion[_,_,{___,V,___},{i_,j_}]/; j < idx
		}];
		result = DeleteCases[result,possibilities],
		{i,Length[lt]}];
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
	A = Prod@@amb[[2]];
	C = Prod@@amb[[3]];
	If[amb[[0]]=== Overlap,
			(*Overlap[ABC,C,A]*)
			{Prod[fi[[2]],C] - Prod[A,fj[[2]]],
				{{A,ToPoly[fj],Prod[]},{-Prod[],ToPoly[fi],C}}},
			(*Inclusion[CBA,C,A]*)
			{fi[[2]] - Prod[A,fj[[2]],C],
			{{A,ToPoly[fj],C},{-Prod[],ToPoly[fi],Prod[]}}}
	]
]


(* ::Text:: *)
(*Same as Poly but without returning the linear combination.*)


SPoly2[amb:_Overlap|_Inclusion,fi_,fj_]:=
	If[amb[[0]]=== Overlap,
			(*Overlap[ABC,C,A]*)
			Prod[fi[[2]],amb[[3]]] - Prod[amb[[2]],fj[[2]]],
			(*Inclusion[CBA,C,A]*)
			fi[[2]] - Prod[amb[[2]],fj[[2]],amb[[3]]]
	]


(* ::Subsection::Closed:: *)
(*Groebner basis*)


(* ::Text:: *)
(*Implementation of the Buchberger algorithm to compute a (partial) Groebner basis of the ideal 'ideal' with at most 'maxiter' iterations being executed (default: 10). The return value is a  set of polynomials {f1,...,fn,g1,...gm} consisting of the elements f1,....,fn from 'ideal' and new elements g1,...,gm. For every new element g, a pair {g, l} is saved in the list cofactors, where l is a list forming a linear combination of g consisting of elements from 'ideal' and certain cofactors.*)
(**)
(*OptionPattern:*)
(*	- Criterion (default: True): Tries to detect and delete redundant ambiguities during the Groebner basis computation.*)
(*	- Ignore (default: 0): A non-negative integer that determines how many elements of the input will be ignored during the first computation of the ambiguities. *)
(*	- MaxDeg (default: Infinity): Only ambiguities with degree smaller than or equal to MaxDeg will be considered during the Groebner basis computation (larger ambiguities are simply ignored). *)
(*	- Info (default: True): Prints information about the computation progress.*)
(*	- Parallel (default: True): Determines whether the computations for which it is possible, are executed in parallel (which speeds up the computation) or in series.*)
(*	- Sorted (default: True):  Sorts the ambiguities before processing in ascending order. This speeds up the computation but results in a different (partial) Groebner basis.*)
(*	- OutputProd (default: False): If this OptionPattern is set to True, the output, i.e. the Groebner basis and the list of cofactors, is given in the Prod data structure. Otherwise, Mathematica's*)
(*	non-commutative multiplication is used.*)
(*	- Rewrite (default: True): Determines whether the cofactors are rewritten in terms of the generators of the ideal. If not, the cofactors consist of all elements of the returned Groebner basis.*)
(**)


SetAttributes[Groebner,HoldFirst]

Groebner[cofactors_,ideal_, maxiter:_?IntegerQ:10, OptionsPattern[{Criterion->True,Ignore->0,MaxDeg->Infinity,Info->False,Parallel->True,Sorted->True,OutputProd->False,Rewrite->True,IterCount->0}]]:=
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
		If[info,Print["Iteration ",count + OptionValue[IterCount], " finished. G has now ", Length[G]," elements\n"]];
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
(*CheckResolvability returns all S-polynomials from the reduction system sys which can not be reduced to zero. Additionally, for each S-polynomial a list containing the linear combination how the S-polynomial was generated from the elements of sys is returned.*)
(*For a description of the OptionPatterns see the documentation of the Groebner method.*)


CheckResolvability[sys_,oldlength:_?IntegerQ:0,OptionsPattern[{Criterion->True,MaxDeg->Infinity,Info->False,Parallel->True,Sorted->True}]]:=
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
		amb = DeleteRedundant[amb,words[[All,1]],Info->info];
	];
	If[OptionValue[Sorted],amb = SortBy[amb,Length[#[[1]]]&]];
	
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
(*F4*)


(* ::Text:: *)
(*Implementation of the Faugere's F4 algorithm to compute a (partial) Groebner basis of the ideal 'ideal' with at most 'maxiter' iterations being executed (default: 10). The return value is a  set of polynomials {f1,...,fn,g1,...gm} consisting of the elements f1,....,fn from 'ideal' and new elements g1,...,gm. For every new element g, a pair {g, l} is saved in the list cofactors, where l is a list forming a linear combination of g consisting of elements from 'ideal' and certain cofactors.*)
(**)
(*OptionPattern:*)
(*	- N (default:100): The number of critical polynomials which are reduced in one step.*)
(*	- Criterion (default: True): Tries to detect and delete redundant ambiguities during the Groebner basis computation.*)
(*	- Ignore (default: 0): A non-negative integer that determines how many elements of the input will be ignored during the first computation of the ambiguities. *)
(*	- MaxDeg (default: Infinity): Only ambiguities with degree smaller than or equal to MaxDeg will be considered during the Groebner basis computation (larger ambiguities are simply ignored). *)
(*	- Info (default: True): Prints information about the computation progress.*)
(*	- Parallel (default: True): Determines whether the computations for which it is possible, are executed in parallel (which speeds up the computation) or in series.*)
(*	- Sorted (default: True):  Sorts the ambiguities before processing in ascending order. This speeds up the computation but results in a different (partial) Groebner basis.*)
(*	- OutputProd (default: False): If this OptionPattern is set to True, the output, i.e. the Groebner basis and the list of cofactors, is given in the Prod data structure. Otherwise, Mathematica's*)
(*	non-commutative multiplication is used.*)
(*	- Rewrite (default: True): Determines whether the cofactors are rewritten in terms of the generators of the ideal. If not, the cofactors consist of all elements of the returned Groebner basis.*)
(**)


SetAttributes[F4,HoldFirst];

F4[cofactors_,ideal_, maxiter:_?IntegerQ:10, OptionsPattern[{N->100,Criterion->True,Ignore->0,MaxDeg->Infinity,Info->False,Parallel->True,Sorted->True,OutputProd->False,Rewrite->True}]]:=
Module[{count,spol,lt,info,G,t1,t2,sorted,oldlength,parallel,maxdeg,n,L,cofactorsL,lc,a,b,rules},
	info = OptionValue[Info];
	sorted = OptionValue[Sorted];
	parallel = OptionValue[Parallel];
	maxdeg = OptionValue[MaxDeg];
	
	cofactors = {};

	G = ToProd/@ideal;
	G = DeleteCases[G,0];
	MakeMonic[G];

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
		RewriteGroebner[cofactors,Info->info];
	][[1]];
	If[info, Print["Rewriting the cofactors took in total ", t2]];
	ToNonCommutativeMultiply[G]
]


(* ::Subsubsection:: *)
(*Symbolic Preprocessing & Reduction*)


SetAttributes[SymbolicPreprocessing,HoldFirst]

SymbolicPreprocessing[{cofactorsF_,columns_},L_,G_]:=
Module[{F,T,lt,g,rules,a,b,newTerms},
	F = L;
	T = Complement[DeleteDuplicates[Flatten[Monomials/@F,1]],columns];
	lt = (LeadingTerm/@G)[[All,2]];
	rules = Table[With[{x = i},{a___,Sequence@@lt[[i]],b___}:>(AppendTo[cofactorsF,{Prod[a],G[[x]],Prod[b]}];Prod[a,G[[x]],b])],{i,Length[G]}];
	While[Length[T] > 0,
		columns = Join[columns,T];
		g = DeleteCases[T/.rules,_List];
		F = Join[F,g];
		newTerms = DeleteDuplicates[Flatten[Monomials/@g,1]];
		T = Complement[newTerms,columns];
	];
	lt = MakeMonic[F];
	cofactorsF = MapIndexed[{#1[[1]]/lt[[#2[[1]]]],#1[[2]],#1[[3]]}&,cofactorsF];
	cofactorsF = DeleteDuplicatesBy[cofactorsF,Prod[#]&];
	DeleteDuplicates[F]
]


Reduction[L_,G_]:=
Module[{F,M,lt,columns,FPlus,a,c,t1,t2,t3,t4,cofactorsF,A,cofactors,pos,f,rule},
	t1 = AbsoluteTiming[
	cofactorsF = L[[All,2]];
	columns = DeleteDuplicates[(LeadingTerm/@L[[All,1]])[[All,2]]];
	F = SymbolicPreprocessing[{cofactorsF,columns},L[[All,1]],G];
	][[1]];
	
	t2 = AbsoluteTiming[
		lt = (LeadingTerm/@F)[[All,2]];
		(*sort in descending order*)
		columns = Reverse[Sort[columns,SortedQ]];
		M = SparseArray[Flatten[MapIndexed[#1/.{Plus->List,c_*Prod[a___]->({#2[[1]],{a}}->c),Prod[a___]->({#2[[1]],{a}}->1)}&,F]/.MapIndexed[#1->#2[[1]]&,columns]]];
	][[1]];
	
	t3 = AbsoluteTiming[
		{A,M} = HermiteDecomposition[M];
		FPlus = DeleteCases[M.(ToProd/@columns),0];
	][[1]];
	
	t4 = AbsoluteTiming[
		cofactors = Map[{}&,cofactorsF];
		Do[
			f = cofactorsF[[rule[[1,2]]]];
			AppendTo[cofactors[[rule[[1,1]]]],{rule[[2]]*f[[1]],f[[2]],f[[3]]}];
			,{rule,ArrayRules[A][[;;-2]]}
		];
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
			(*Overlap[ABC,A,C]*)
			If[Prod[fi,amb[[3]]] === Prod[amb[[2]],fj], 
				Sequence@@{},
				Sequence@@{{Prod[fi,amb[[3]]],{Prod[],fi,Prod[amb[[3]]]}},{Prod[amb[[2]],fj],{Prod[amb[[2]]],fj,Prod[]}}}
			],
			(*Inclusion[ABC,A,C]*)
			If[fi === Prod[amb[[2]],fj,amb[[3]]], 
				Sequence@@{},
				Sequence@@{{fi,{Prod[],fi,Prod[]}},{Prod[amb[[2]],fj,amb[[3]]],{Prod[amb[[2]]],fj,Prod[amb[[3]]]}}}
			]
	]


CheckResolvabilityF4[G_,oldlength:_?IntegerQ:0,OptionsPattern[{Criterion->True,MaxDeg->Infinity,Info->False,Parallel->True,Sorted->True}]]:=
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
		amb = DeleteRedundant[amb,words[[All,1]],Info->info]
	];
	If[sorted,amb = SortBy[amb,Length[#[[1]]]&]];
	
	(*generate S-polynomials*)
	t2 = AbsoluteTiming[
		If[parallel,
			spol = ParallelMap[SPolyF4[#,G[[#[[4,1]]]],G[[#[[4,2]]]]]&,amb,DistributedContexts->Automatic,Method->"CoarsestGrained"],
			spol = Map[SPolyF4[#,G[[#[[4,1]]]],G[[#[[4,2]]]]]&,amb]
		];		
	][[1]];
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
(*Groebner basis without cofactors*)


(* ::Text:: *)
(*Implementation of the Buchberger algorithm to compute a (partial) Groebner basis of the ideal 'ideal' with at most 'maxiter' iterations being executed (default: 10). The return value is a  set of polynomials {f1,...,fn,g1,...gm} consisting of the elements f1,....,fn from 'ideal' and new elements g1,...,gm. *)
(**)
(*OptionPattern:*)
(*	- Criterion (default: True): Tries to detect and delete redundant ambiguities during the Groebner basis computation.*)
(*	- Ignore (default: 0): A non-negative integer that determines how many elements of the input will be ignored during the first computation of the ambiguities. *)
(*	- MaxDeg (default: Infinity): Only ambiguities with degree smaller than or equal to MaxDeg will be considered during the Groebner basis computation (larger ambiguities are simply ignored). *)
(*	- Info (default: True): Prints information about the computation progress.*)
(*	- Parallel (default: True): Determines whether the computations for which it is possible, are executed in parallel (which speeds up the computation) or in series.*)
(*	- Sorted (default: True):  Sorts the ambiguities before processing in ascending order. This speeds up the computation but results in a different (partial) Groebner basis.*)
(*	*)


GroebnerWithoutCofactors[ideal_,maxiter:_?IntegerQ:10,OptionsPattern[{Ignore->0, MaxDeg->Infinity,Info->False,Parallel->True,Sorted->True,Criterion->True}]]:=
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


CheckResolvability2[sys_,oldlength:_?IntegerQ:0,OptionsPattern[{Criterion->True,MaxDeg->Infinity,Info->False,Sorted->True,Parallel->True}]]:=
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
		amb = DeleteRedundant[amb,words[[All,1]],Info->info]
	];
	If[OptionValue[Sorted],amb = SortBy[amb,Length[#[[1]]]&]];
	
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
	DeleteCases[Map[CreateRedSys,polies],{}]


CreateRedSys[p_]:=
Module[{lt,poly},
	poly = ToProd[p];
	If[poly===0,
		Return[{}]
	];
	lt = LeadingTerm[poly];
	{lt[[2]], Expand[-1/lt[[1]]*Remainder[poly,lt]]}
]


ToPoly[poly:{_List,_Prod|_Plus|_Times}]:= Prod@@poly[[1]] - poly[[2]]


ToPoly[poly:{_List,0}]:= Prod@@poly[[1]]


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


IsLinearCombination[triples_List, G_List] := AllTrue[triples[[All,2]],MemberQ[G,#]&]


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
	While[i <= s,
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
(*Adjungate operator definition & auxiliary stuff*)


(* ::Text:: *)
(*Properties of the adjoint operator. Compatible with Prod[] and Mathematica's non-commutative multiplication.*)


adj[Prod[a_]]:=Prod[adj[a]]
adj[Prod[a__,b__]] := Prod[adj[Prod[b]],adj[Prod[a]]]
adj[Prod[]]:= Prod[]
adj[Times[a__,Prod[b___]]]:= a adj[Prod[b]]
adj[a___,b_Plus,c___]:=(adj[a,#,c]&/@b)
adj[adj[a__]] := a


adj[NonCommutativeMultiply[a_,b_]] := NonCommutativeMultiply[adj[b],adj[a]]
adj[a_?CoeffQ]:= a
adj[-a_]:= -adj[a];
adj[Times[a_?CoeffQ,NonCommutativeMultiply[b___]]]:= a adj[NonCommutativeMultiply[b]]
adj[Times[a_?CoeffQ,b_]]:=a adj[b]


Pinv[a_] := {a**SuperDagger[a]**a-a,SuperDagger[a]**a**SuperDagger[a]- SuperDagger[a],adj[SuperDagger[a]]**adj[a]-a**SuperDagger[a],adj[a]**adj[SuperDagger[a]]-SuperDagger[a]**a};
Pinv[a_,b_] := {a**b**a-a,b**a**b- b,adj[b]**adj[a]-a**b,adj[a]**adj[b]-b**a};


AddAdj[S_List] := Join[S,adj/@S];


(* ::Subsection::Closed:: *)
(*Quiver*)


(* ::Text:: *)
(*Data structure*)


Quiver = {RepeatedNull[{__,__,__}]};


(* ::Text:: *)
(*Gives all Sources and Targets, respectively, of a certain label of a quiver Q (not necessarily with unique labels).*)


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
(*Returns the signature of a polynomial w.r.t. to a quiver Q (not necessarily with unique labels)*)


QSignature[l_List,Q:Quiver]:= Map[QSignature[#,Q]&,l]


QSignature[p_,Q:Quiver]:=
Module[{monomials,m,results,i,begin,end,comb,vertices,x},
	monomials = (MonomialList[ToProd[p]]/.c__*Prod[x___]->Prod[x])/.Prod->List;
	Sort[Catch[
	(*base case: only one mononmial*)
	If[Length[monomials] === 1,
			m = monomials[[1]];
			(*zero has signature V x V*)
			If[m===0,
				vertices = Sort[DeleteDuplicates[Flatten[Q[[All,2;;3]]]]];
				Throw[Tuples[vertices,2]]
			];
			(*empty monomial = Prod[] \[Rule] return all empty paths*)
			If[Length[m]===0, 
				vertices = Sort[DeleteDuplicates[Flatten[Q[[All,2;;3]]]]];
				Throw[Map[{#,#}&,vertices]]
			];
			(*get all sources and targets of first label*)
			begin = Cases[Q,{m[[-1]],_,_}][[All,2;;3]];
			For[i = Length[m]-1, i > 0, i--,
				(*get all sources and targets of label i*)
				end = Cases[Q,{m[[i]],_,_}][[All,2;;3]];
				(*check if the two labels are compatible*)
				comb = {};
				Do[Map[If[i[[2]]===#[[1]],AppendTo[comb,{i[[1]],#[[2]]}]]&,end],{i,begin}];
				If[comb === {}, Throw[{}]];
				begin = DeleteDuplicates[comb];
			];
			Throw[begin],	
			(*usual case: split polynomial into monomials*)
			results = Map[QSignature[Prod@@#,Q]&,monomials];
			Throw[Intersection[Sequence@@results]]
	]
	]]
]


CompatibleQ[p_,Q:Quiver]:= QSignature[p,Q]=!={}


UniformlyCompatibleQ[p_,Q:Quiver]:= 
Module[{monomials,signatures},
	If[!CompatibleQ[p,Q],
		Return[False]
	];
	monomials = MonomialList[p];
	signatures = Map[QSignature[#,Q]&,monomials];
	AllTrue[signatures,#===signatures[[1]]&]
]


PlotQuiver[Q:Quiver]:=
	GraphPlot[Map[{#[[2]]->#[[3]],#[[1]]}&,Q],DirectedEdges->True,SelfLoopStyle->.2]


TrivialQuiver[F_]:= Module[{vars,v},
	vars = DeleteDuplicates[Select[(ToProd/@F)/.{Prod->List,Plus->List,Times->List}//Flatten,!CoeffQ[#]&]];
	Table[{v,1,1},{v,vars}]
]


QOrderCompatibleQ[p_,Q_]:= Module[{lm},
	lm = ToProd[LeadingTerm[p][[2]]];
	QSignature[lm,Q] === QSignature[p,Q]
]


QConsequenceQ[certificate_,F_,Q_]:= Module[{f,signaturef},
	(* check if certificate is indeed a liner combination of elements in F *)
	If[!IsLinearCombination[certificate,F],
		Print["The input is not a linear combination of elements of the given set."];
		Return[False]
	];  

	f = MultiplyOut[certificate];
	signaturef = QSignature[f,Q];
	
	(* f has to be compatible *)
	If[signaturef === {}, Return[False]];
	
	(* Q-consequence test has to pass for each summand in the certificate *)
	!MemberQ[Map[QConsequenceTest[#,signaturef,Q]&,certificate],False]
]


QConsequenceTest[{ai_,gi_,bi_},signaturef_,Q_]:= Module[{sigGi,sigMonomials,mi},
	(* check if there is m_i \in \supp(g_i) s.t. \sigma(m_i) = \sigma(g_i) *)
	sigGi = QSignature[gi,Q];
	sigMonomials = QSignature[MonomialList[gi],Q];
	If[!MemberQ[sigMonomials,sigGi], Return[False]];
	
	(* check if \sigma(f) \subseteq \sigma(a_i m_i b_i) *)
	mi = Extract[MonomialList[gi],FirstPosition[sigMonomials,sigGi,1]];
	SubsetQ[QSignature[Prod[ai,mi,bi],Q],signaturef]
]


(* ::Subsubsection:: *)
(*Q-completion*)


SetAttributes[QCompletion,HoldFirst]

QCompletion[cofactors_,ideal_, maxiter:_?IntegerQ:10, OptionsPattern[{Criterion->True,Ignore->0,MaxDeg->Infinity,Info->False,Parallel->True,Sorted->True,OutputProd->False,Rewrite->True,IterCount->0}]]:=
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
		If[info,Print["Iteration ",count + OptionValue[IterCount], " finished. G has now ", Length[G]," elements\n"]];
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


(* ::Subsection::Closed:: *)
(*Certify*)


Certify[assumptionsInput_List,claimsInput_,Q:Quiver,OptionsPattern[{MaxIter->10,MaxDeg->Infinity,MultiLex->False,Info->False,Parallel->True,Sorted->True,Criterion->True}]]:=
 Module[{info,maxiter,reduced,vars,cofactors,G,sigAssump,sigClaim,certificate,rules,lc,toIgnore,toIgnoreOld,zeros,i,knowns,unknowns,t,assumptions,claims,redCofactors,k,l,count,assumptionsRed,a,b},
	info = OptionValue[Info];
	maxiter = OptionValue[MaxIter];
	
	assumptions = ToProd/@assumptionsInput;
	assumptions = DeleteCases[assumptions,0];
	
	If[Head[claimsInput] === List,
		claims = claimsInput,
		claims = {claimsInput}
	];

	(*check compatibility of the assumptions and the claims*)
	sigAssump = Map[UniformlyCompatibleQ[#,Q]&,assumptions];
	If[MemberQ[sigAssump,False],
		Print["The assumption ", Extract[assumptionsInput,Position[sigAssump,False][[1,1]]]," is not uniformly compatible with the quiver."]; Return[$Failed]
	];
	sigClaim = Map[QSignature[#,Q]&,claims];
	If[MemberQ[sigClaim,{}], 
		Print["The claim ", Extract[claims,Position[sigClaim,{}][[1,1]]]," is not compatible with the quiver."]; Return[$Failed]	
	];

	(*set up the ring*)
	If[info,
		Print["Using the following monomial ordering:"]
	];
	If[OptionValue[MultiLex],
		knowns = DeleteDuplicates[Cases[Q[[All,1]],Alternatives@@Flatten[Map[#/.Alternatives[Plus,Times,NonCommutativeMultiply]->List&,claims]]]];
		unknowns = DeleteDuplicates[Cases[Q[[All,1]],var_/;!MemberQ[knowns,var]]];
		SetUpRing[knowns,unknowns,Info->info],
		SetUpRing[DeleteDuplicates[Q[[All,1]]],Info->info]
	];

	(*make ideal monic*)
	lc = MakeMonic[assumptions];
	
	(*interreduce the generators*)
	{assumptionsRed,redCofactors} = Interreduce[assumptions,InputProd->True];
	If[info, Print["\nInterreduced the input from ", Length[assumptions], " polynomials to ", Length[assumptionsRed], ".\n"]];
	
	(*compute the Groebner basis and reduce the claims*)
	If[info, Print["Computing a (partial) Groebner basis and reducing the claim...\n"]];
	(*do computation iteratively*)
	cofactors = {};
	zeros = ConstantArray[0,Length[claims]];
	toIgnore = Length[assumptionsRed];
	i = 1;
	If[info,Print["Starting iteration ", i ,"...\n"]];
	G = Groebner[cofactors,assumptionsRed,1,MaxDeg->OptionValue[MaxDeg],Info->OptionValue[Info],Parallel->OptionValue[Parallel],Sorted->OptionValue[Sorted],Criterion->OptionValue[Criterion],OutputProd->True,Rewrite->False,IterCount->i-1];
	i++;
	reduced = ReducedForm[vars,G,claims];
	While[reduced =!= zeros && i <= maxiter,
		toIgnoreOld = Length[G];
		If[info,Print["Starting iteration ", i ,"...\n"]];
		G = Groebner[cofactors,G,1,Ignore->toIgnore,MaxDeg->OptionValue[MaxDeg],Info->OptionValue[Info],Parallel->OptionValue[Parallel],Sorted->OptionValue[Sorted],Criterion->OptionValue[Criterion],OutputProd->True,Rewrite->False,IterCount->i-1];
		i++;
		toIgnore = toIgnoreOld;
		count = 0;
		While[reduced =!= zeros && count < 4,
			reduced = ReducedForm[vars,RandomSample[G],claims];
			count++;
		];
	];
	
	(*rewrite the linear combination*)
	If[info, Print["Rewriting the cofactors has started..."]];
	t = AbsoluteTiming[certificate = RewriteCertify[vars,cofactors];][[1]];
	If[info, Print["Rewriting the cofactors took in total ", t]];
	
	(*rewrite in terms of assumptions and not of the interreduced assumptions*)
	rules = Table[{k_,assumptionsRed[[i]],l_}->
		Sequence@@Table[{Prod[k,j[[1]]],j[[2]],Prod[j[[3]],l]},{j,redCofactors[[i]]}],{i,Length[redCofactors]}];
	certificate = certificate/.rules;
	(*take care of leading coefficients in the certificate*)
	rules = MapIndexed[{a_,#1,b_}->{a/lc[[First[#2]]],Expand[lc[[First[#2]]]*#1],b}&,assumptions];	
	certificate = certificate/.rules;
	(*convert back to NonCommutativeMultiply*)
	
	certificate = ToNonCommutativeMultiply[Map[ExpandLeft[CollectLeft[#]]&,certificate]];
	
	(*return the reduced claims and the linear combinations*)
	If[Head[claimsInput] =!= List, 
		sigClaim = sigClaim[[1]]; reduced = reduced[[1]]; certificate = certificate[[1]]; zeros = zeros[[1]];
	];
	
	If[info,
		(*full info*)
		If[reduced === zeros,
			Print["\nDone! All claims were successfully reduced to 0."],
			Print["\nDone! Not all claims could be reduced to 0."]
		];
		{Map[QSignature[#,Q]&,assumptions],sigClaim,reduced,certificate},
		(*only basic info*)
		If[reduced === zeros, certificate, $Failed]
	]
]


RewriteCertify[varsInput_,cofactors_]:= Module[
{a,b,i,j,count,rules,toReduce,occurring,toAdd,g,vars},
	If[MatchQ[varsInput,{{{_,_,_}...}...}],
		toReduce = Map[ToProd,varsInput,{3}];
		vars = Flatten[toReduce,1],
		toReduce = Map[ToProd,varsInput,{2}];
		vars = toReduce
	];
				
	toAdd = Position[cofactors[[All,1]],Alternatives@@vars[[All,2]],{1}];
	occurring = {};
	(*find all cofactors actually appearing in the certificate*)
	While[toAdd =!= {},
		occurring = Join[occurring,toAdd];
		g = Extract[cofactors,toAdd];
		toAdd = Flatten[Map[Position[cofactors[[All,1]],Alternatives@@#[[2,All,2]],{1}]&,g],1];
		toAdd = Complement[toAdd,occurring];
	];
	If[occurring === {},
		Return[toReduce]
	];
	occurring = Sort[occurring];
	
	(*sucessively make reduction rules*)
	g = cofactors[[occurring[[1,1]]]];
	rules = {{a__,g[[1]],b__} -> Sequence@@Table[{Prod[Evaluate[a],j[[1]]],j[[2]],Prod[j[[3]],Evaluate[b]]},{j,g[[2]]}]};
	Do[
		g = f[[2]]/.rules;
		AppendTo[rules,{a__,f[[1]],b__} -> Sequence@@Table[{Prod[Evaluate[a],j[[1]]],j[[2]],Prod[j[[3]],Evaluate[b]]},{j,g}]];
	,{f,Extract[cofactors,occurring[[2;;]]]}];
	
	toReduce/.rules
]


(* ::Subsection::Closed:: *)
(*End*)


Copyright[a_String,b___String]:=Print[StringJoin[Prepend[{"\n",#}&/@{b},a]]]


Copyright[
    "Package OperatorGB version 1.1.1",
    "Copyright 2019, Institute for Algebra, JKU",
    "written by Clemens Hofstadler"];


End[]


EndPackage[]


?OperatorGB
