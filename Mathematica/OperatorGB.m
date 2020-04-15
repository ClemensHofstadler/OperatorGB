(* ::Package:: *)

BeginPackage["OperatorGB`"]


Global`OperatorGB::usage="";


Clear[
	CoeffQ,Prod,ToProd,ToNonCommutativeMultiply,
	SetUpRing,
	WordOrder,varSets,
	LeadingTerm,DegLex,WeightedDegLex,MultiLex,Weight,SortedQ,
	ExtractReducibleWords,GenerateAmbiguities,Overlap,Inclusion,DeleteRedundant,
	Groebner,CheckResolvability,
	F4,
	ReducedForm,
	GroebnerWithoutCofactors,ApplyRules,
	CreateRedSys,ToPoly,Rewrite,Interreduce,
	MultiplyOut,LinearCombinationQ,CertificateCoeffQ,IntegerCoeffQ,CheckCertificate,CheckCertificates,
	CollectLeft,CollectRight,ExpandLeft,ExpandRight,
	adj,Pinv,AddAdj,IntegerCoeffQ,
	Quiver,QSignature,PlotQuiver,CompatibleQ,UniformlyCompatibleQ,QOrderCompatibleQ,QConsequenceQ,QConsequenceQCrit,
	QCompletion,
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
Overlap::usage="Datastructure for overlap ambiguities."
Inclusion::usage="Datastructure for inclusion ambiguities."
ExtractReducibleWords::usage="Preprocessing before GenerateAmbiguities."
GenerateAmbiguities::usage="GenerateAmbiguities[words] computes all ambigutites among all words in the set 'words'."
DeleteRedundant::usage="Chain criterion to remove redundant ambiguities."


(*Groebner basis*)
Groebner::usage="Groebner[cofactors_,ideal_, maxiter:_?IntegerQ:10, OptionsPattern[{Criterion->True,Ignore->0,MaxDeg->Infinity,Info\[Rule]False,Parallel->True,Sorted->True}]] executes at most maxiter iterations of the Buchberger algorithm to compute
a (partial) Groebner basis of an ideal. Additionally, for every new element in the Groebner basis a list of cofactors is saved in the list cofactors forming a linear combination of the new element. For further information concerning the OptionPatterns
please see the documentation or the source code."
CheckResolvability::usage=""


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
Rewrite::usage="Rewrite[cofactorsF, cofactorsG] rewrites a linear combination, which is safed in cofactorsF, with the elements from cofactorsG."
Interreduce::usage="Interreduce[ideal_] interreduces the polynomials in 'ideal'."


(*Check certificate*)
MultiplyOut::usage="To multiply out a list of cofactors given in terms of the built in non-commutative multiplication."
LinearCombinationQ::usage="Checks whether a given set of triples is a linear combination of a set of polynomials."
CertificateCoeffQ::usage=""
IntegerCoeffQ::usage="Checks whether a given certificate only contains integer coefficients"
CheckCertificate::usage="Check that a single certificate indeed gives the claim, is w.r.t. to the assumptions and that it only consists of integer coefficients"
CheckCertificates::usage="Applies CheckCertificate to several certificates"


(*Collect cofactors*)
CollectLeft::usage="Tries to collect triples in a cofactor representation having the same left cofactors."
CollectRight::usage="Tries to collect triples in a cofactor representation having the same right cofactors."
ExpandLeft::usage=""
ExpandRight::usage=""


(*Adjungate operator definition & auxiliary stuff*)
adj::usage="adj[A] represents the adjoint of the operator A"
Pinv::usage="Gives the 4 Moore-Penrose equations."
AddAdj::usage="Adds the adjoint statements to a given list of statements."
IntegerCoeffQ::usage="Tests if a certficate contains only integer coefficients."


(*Quiver*)
Quiver::usage="Data structure of a Quiver"
QSignature::usage="QSignature[poly,Q] returns the signature of the polynomial poly w.r.t. the quiver Q (not necessarily with unique labels)"
PlotQuiver::usage="Plot a quiver Q"
CompatibleQ::usage="Tests whether a polynomial is compatible with a quiver."
UniformlyCompatibleQ::usage="Tests whether a polynomial is uniformly compatible with a quiver."
TrivialQuiver::usage="Returns the trivial quiver containing all variables of the given input."
QOrderCompatibleQ::usage="Tests whether a polynomial is Q-order-compatible with a quiver and the order defined by SetUpRing."
QConsequenceQ::usage="Tests whether a certificate is a Q-consequence of a set of polynomials and a quiver using the definition of Q-consequence."
QConsequenceQCrit::usage="Tests whether a certificate is a Q-consequence of a set of polynomials and a quiver using a criterion."


(*Q-completion*)
QCompletion::usage="Q-completion procedure"


(*Certify*)
Certify::usage="Certifies whether a certain claim is a consequence of some assumptions via Groebner 
basis computations. Additionally, compatibility with a given quiver is checked."


Begin["`Private`"]


(* ::Subsection::Closed:: *)
(*Noncommutative multiplication*)


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


LeadingTerm[poly_] := LeadingTermIntern[ToProd[poly]]


LeadingTermIntern[poly_Prod]:=
	List[1,List@@poly]


LeadingTermIntern[poly_]:=
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


GenerateAmbiguities[l_List,newpart_List,maxdeg_:Infinity,OptionsPattern[Parallel->True]]:= 
	Select[Join[GenerateOverlaps[l,newpart,Parallel->OptionValue[Parallel]],
		GenerateInclusions[l,newpart,Parallel->OptionValue[Parallel]],
		GenerateAmbiguities[newpart,maxdeg,Parallel->OptionValue[Parallel]]], Length[#[[1]]] <= maxdeg &]


GenerateAmbiguities[l_List,maxdeg_:Infinity,OptionsPattern[Parallel->True]] := 
	Select[Join[GenerateOverlaps[l,Parallel->OptionValue[Parallel]],
		GenerateInclusions[l,Parallel->OptionValue[Parallel]]], Length[#[[1]]] <= maxdeg &]


(* ::Text:: *)
(*Chain Criterion*)


DeleteRedundantSimple[amb_List,lt_List,OptionsPattern[{Info->False}]]:= 
Module[{result,pattern,t,i,j,k,V,l},
	t = AbsoluteTiming[
	result = amb;
	Do[
		V = Sequence@@l[[1]];
		k = l[[2]];
		pattern = With[{idx = k}, _[{___,V,___},_,_,{i_,j_}]/; idx < i && idx < j];
		result = DeleteCases[result,pattern]
	,{l,lt}];
	][[1]];
	If[OptionValue[Info],
		Print["Removed ", Length[amb] - Length[result], " ambiguities in ",t]];
	result
]


DeleteRedundantComplex[amb_List,lt_List,OptionsPattern[{Info->False}]]:= 
Module[{result,overlaps,incls,t,i,j,pattern,V,A,C,X,Y,pre,post,l,ltsort,k,L,R},
	t = AbsoluteTiming[
	overlaps = Cases[amb,_Overlap];
	incls = Cases[amb,_Inclusion];
	ltsort = SortBy[lt,Length[#[[1]]]&];
	Do[
		V = Sequence@@l[[1]];
		k = l[[2]];
		pattern =  Alternatives@@{
		(* V | A *)
		_[_,{___,V,___},_,_],
		(* V | C *)
		_[_,_,{___,V,___},_],
		(* ABC = LVR and |L| = 0 and |R| < |C| *)
		_[{V,R___},_,C_,_]/; Length[{R}] < Length[C],
		(* ABC = LVR and 0 < |L| < |A| *)
		_[{L__,V,___},A_,_,_]/; Length[{L}] < Length[A],
		(* ABC = LVR and |L| \[GreaterEqual] |A| and |R| > 0 *)
		_[{L__,V,__},A_,_,_]/; Length[{L}] >= Length[A],
		(* V |\[NonBreakingSpace]B *)
		_[{L__,V,R__},A_,C_,_]/;Length[{L}] >= Length[A] && Length[{R}] >= Length[C]
		};
		overlaps = DeleteCases[overlaps,pattern];
	,{l,DeleteDuplicates[ltsort,Length[SequencePosition[#1[[1]],#2[[1]]]] > 0&]}];
	Do[
		V = Sequence@@l[[1]];
		k = l[[2]];
		pattern = {
		(* V | A *)
		_[_,{___,V,___},_,{_,j_}]/; k < j,
		(* V | C *)
		_[_,_,{___,V,___},{_,j_}]/; k < j,
		(* V | B and |AC| > 0 *)
		_[{L__,V,R__},A_,C_,{_,j_}]/; k < j && Length[{L}] >= Length[A] && Length[{R}] >= Length[C] && Length[A] + Length[C] > 0,
		(* B | V and |V| < |ABC| *)
		_[{L__,V,R__},A_,C_,{_,j_}]/; k < j && Length[{L}] < Length[A] && Length[{R}] < Length[C]
		};
		incls = DeleteCases[incls,pattern],
	{l,ltsort}];
	][[1]];
	result = Join[overlaps,incls];
	If[OptionValue[Info],
		Print["Removed ", Length[amb] - Length[result], " ambiguities in ",t]];
	result
]


DeleteRedundant[amb_List,lt_List,OptionsPattern[{Info->False}]]:= 
Module[{result,overlaps,incls,t,i,j,pattern,V,A,C,X,Y,l,k},
	t = AbsoluteTiming[
	overlaps = Cases[amb,_Overlap];
	incls = Cases[amb,_Inclusion];
	Do[
		V = Sequence@@l[[1]];
		k = l[[2]];
		pattern =  _[{___,V,___},_,_,{i_,j_}]/; k < i && k < j;
		overlaps = DeleteCases[overlaps,pattern];
		pattern = Alternatives@@{
		(* V | A *)
		_[_,{___,V,___},_,{_,j_}]/; k < j,
		(* V | C *)
		_[_,_,{___,V,___},{_,j_}]/; k < j,
		(* V | B and |AC| > 0 *)
		_[{X__,V,Y__},A_,C_,{_,j_}]/; k < j && Length[{X}] >= Length[A] && Length[{Y}] >= Length[C] && Length[A] + Length[C]\[NonBreakingSpace]> 0,
		(* B | V | ABC and V \[NotEqual] ABC *)
		_[{X__,V,Y__},A_,C_,{_,j_}]/; k < j && Length[{X}] < Length[A]\[NonBreakingSpace]&& Length[{Y}] < Length[C]
		};
		incls = DeleteCases[incls,pattern],
		{l,lt}];
	][[1]];
	result = Join[overlaps,incls];
	If[OptionValue[Info],
		Print["Removed ", Length[amb] - Length[result], " ambiguities in ",t]];
	result
]


(* ::Text:: *)
(*Gebauer-M\[ODoubleDot]ller criterion*)


GebauerMoeller[ambInput_List,lt_List,OptionsPattern[{Info->False}]]:= 
Module[{amb,result,t,A,C,i,s,pattern,a,j,APrime,CPrime,pos,selected,idx},
	t = AbsoluteTiming[
	amb = SortBy[ambInput,Length[#[[1]]]&];

	If[Length[amb] === 0,
		result = {},
		idx = Flatten[amb[[All,4]]];
		result = Reap[
			Do[
			selected = Select[amb,Max[#[[4]]] === s&];
			While[Length[selected]\[NonBreakingSpace]> 0,
				a = First[selected];
				selected = Delete[selected,{1}];
				Sow[a];
			
				APrime = Sequence@@a[[2]]; CPrime = Sequence@@a[[3]];
				j = Min[a[[4]]];
				Switch[a,
				Overlap[_,_,_,{s,_}],
					pattern = Alternatives@@{
						Overlap[_,_,{CPrime,___},{s,i_}]/; i > j,
						Overlap[_,_,{CPrime,__},{s,i_}]/; i <= j,
						Overlap[_,A_,{CPrime},{s,i_}]/; i === j && Length[A] > Length[{APrime}],
						Inclusion[_,{},{CPrime,___},{i_,s}]/; i > j,
						Inclusion[_,{},{CPrime,__},{i_,s}]/; i <= j
			    },
				Overlap[_,_,_,{_,s}], 
			    pattern = Alternatives@@{
						Overlap[_,{___,APrime},_,{i_,s}]/; i > j,
						Overlap[_,{__,APrime},_,{i_,s}]/; i <= j,
						Inclusion[_,{___,APrime},{},{i_,s}]/; i > j,
						Inclusion[_,{__,APrime},{},{i_,s}]/; i <= j
			    },
				Inclusion[_,_,_,{s,_}], 
					pattern = Alternatives@@{
						Overlap[_,_,_,{i_,s}]/; i > j,
						Overlap[_,_,{__},{i_,s}]/; i <= j,
						Overlap[_,_,_,{s,_}],
						Inclusion[_,_,_,{i_,s}]/; i > j,
						Inclusion[_,A_,C_,{i_,s}]/; i <= j && Length[A] + Length[C]\[NonBreakingSpace]> 0,
						Inclusion[_,_,_,{s,i_}]/; i > j,
						Inclusion[_,_,_,{s,i_}]/; i === j && Length[A] > Length[{APrime}]
				},
				_,
					pattern = Alternatives@@DeleteCases[{
						If[Length[{CPrime}] > 0,
							Sequence@@{
								Overlap[_,{___,APrime},_,{i_,s}]/; i > j,
								Overlap[_,{__,APrime},_,{i_,s}]/; i <= j
							},
							{}],
						If[Length[{APrime}] > 0,
							Sequence@@{
								Overlap[_,_,{CPrime,___},{i_,s}]/; i >= j,
								Overlap[_,_,{CPrime,__},{i_,s}]/; i < j
							},
							{}],
						Inclusion[_,{___,APrime},{CPrime,___},{i_,s}]/; i > j,
						Inclusion[_,{__,APrime},{CPrime,__},{i_,s}]/; i <= j,
						If[Length[{APrime}] > 0 && Length[{CPrime}] > 0,
							Sequence@@{
								Inclusion[_,_,_,{s,i_}]/; i > j,
								Inclusion[_,{__},_,{s,i_}]/; i <= j
							},
							{}]	
					},{}];
				];
				selected = DeleteCases[selected,pattern];
			],
			{s,Min[idx],Max[idx]}];
			][[2,1]];
	]
	][[1]];
	If[OptionValue[Info],
		Print["Removed ", Length[ambInput] - Length[result], " ambiguities in ",t]];
	Join[Cases[result,_Overlap],Cases[result,_Inclusion]]
]


(* ::Subsection::Closed:: *)
(*S-Polynomials*)


(* ::Text:: *)
(*SPoly[amb,fi,fj] computes the S-polynomial corresponding to the ambiguity amb, which comes from the two polynomials fi and fj. Additionally, the linear combination how the S-polynomial was computed from fi and fj is returned in a list. *)


SPoly[amb:_Overlap|_Inclusion,fi_,fj_]:=
Module[{A,C},
	A = Prod@@amb[[2]];
	C = Prod@@amb[[3]];
	If[amb[[0]]=== Overlap,
			(*Overlap[ABC,C,A]*)
			{Prod[fi[[2]],C] - Prod[A,fj[[2]]],
				{{A,amb[[4,2]],Prod[]},{-Prod[],amb[[4,1]],C}}},
			(*Inclusion[CBA,C,A]*)
			{fi[[2]] - Prod[A,fj[[2]],C],
			{{A,amb[[4,2]],C},{-Prod[],amb[[4,1]],Prod[]}}}
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
(*Implementation of the Buchberger algorithm to compute a (partial) Groebner basis of the ideal 'ideal' with at most 'maxiter' iterations being executed (default: 10). The return value is a  set of polynomials G = {f1,...,fn,g1,...gm} consisting of the elements f1,....,fn from 'ideal' and new elements g1,...,gm. For every element g\in G a list forming a linear combination of g consisting of elements from 'ideal' and certain cofactors is saved in the list cofactors.*)
(**)
(*OptionPattern:*)
(*	- Criterion (default: True): Tries to detect and delete redundant ambiguities during the Groebner basis computation.*)
(*	- Ignore (default: 0): A non-negative integer that determines how many elements of the input will be ignored during the first computation of the ambiguities. *)
(*	- MaxDeg (default: Infinity): Only ambiguities with degree smaller than or equal to MaxDeg will be considered during the Groebner basis computation (larger ambiguities are simply ignored). *)
(*	- Info (default: False): Prints information about the computation progress.*)
(*	- Parallel (default: True): Determines whether the computations for which it is possible, are executed in parallel (which speeds up the computation) or in series.*)
(*	- Sorted (default: True):  Sorts the ambiguities before processing in ascending order. This speeds up the computation but results in a different (partial) Groebner basis.*)
(*	- IterCount (default: 0): defines from which number the iterations are counted (only relevant for the printed information)*)
(**)


SetAttributes[Groebner,HoldFirst]

Groebner[cofactors_,ideal_, maxiter:_?IntegerQ:10, OptionsPattern[{Criterion->True,Ignore->0,MaxDeg->Infinity,Info->False,Parallel->True,Sorted->True,IterCount->0}]]:=
Module[{lc,count,spol,lt,info,p,h,G,r,t1,t2,rules,sorted,oldlength,parallel,hrule,maxdeg,intern,criterion,i},
	info = OptionValue[Info];
	sorted = OptionValue[Sorted];
	parallel = OptionValue[Parallel];
	maxdeg = OptionValue[MaxDeg];
	criterion = OptionValue[Criterion];
	intern = FreeQ[ideal,NonCommutativeMultiply]; 

	If[intern,
		G = ideal,
		G = ToProd/@ideal;
		lc = MakeMonic[G];
		cofactors = MapIndexed[{{1/#1*Prod[],#2[[1]],Prod[]}}&,lc];
	];
	
	G = CreateRedSys[G];
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
				lt = LeadingTermIntern[h];
				If[lt[[1]] =!= 1, 
					h = Expand[1/lt[[1]]*h]; p[[2]] = (ReplacePart[#,1 -> 1/lt[[1]]*#[[1]]]&/@ p[[2]])
				];
				hrule = CreateRedSys[h];
				AppendTo[G,hrule];
				AppendTo[cofactors,p[[2]]]; 
				AppendTo[rules,ExtractRule[hrule,Length[G]]];
			];
			i--;
		,{p,spol}];,i];][[1]];
		
		count++;
		If[info, 
			Print["The second reduction took ", t1];
			Print["Iteration ",count + OptionValue[IterCount], " finished. G has now ", Length[G]," elements\n"]
		];
		If[count < maxiter, 
			spol = CheckResolvability[G,oldlength,Criterion->criterion,MaxDeg->maxdeg,Info->info,Sorted->sorted,Parallel->parallel];
			oldlength = Length[G];
		];
	];
	
	G = ToPoly[G];
	
	If[intern,
		G,
		RewriteGroebner[cofactors,ideal,Info->info];
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
		amb = DeleteRedundant[amb,words,Info->info];
	];
	If[OptionValue[Sorted],
		amb = SortBy[amb,Length[#[[1]]]&];
	];
	
	(*generate S-polynomials*)
	t2 = AbsoluteTiming[
	spol = DeleteCases[
		If[parallel,
			ParallelMap[SPoly[#,sys[[#[[4,1]]]],sys[[#[[4,2]]]]]&,amb,DistributedContexts->Automatic,Method->"CoarsestGrained"],
			Map[SPoly[#,sys[[#[[4,1]]]],sys[[#[[4,2]]]]]&,amb]
		]
	,{0,_}];
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
											])&,spol],{}],#[[1]]&];
	][[1]];
	If[info,Print["Reducing S-polys: ",t3, " (",Length[spol], " remaining)"]];
	If[OptionValue[Sorted],
		If[$VersionNumber < 12,
			Sort[spol,SortedQ[LeadingTermIntern[#1[[1]]][[2]],LeadingTermIntern[#2[[1]]][[2]]]&],
			SortBy[spol,LeadingTermIntern[#[[1]]][[2]]&,SortedQ]
		],
		spol
	]
]


(* ::Text:: *)
(*ExtractRules[vars,sys] generates a set of reduction rules out of the reduction system sys with the special feature that the cofactors of the reduction will be sowed whenever such a rule is applied. They can then be reaped using the Reap command.*)


ExtractRules[sys_]:=
Module[{a,b,coeff,i,p,q,f},
	a=Unique[];b=Unique[];coeff=Unique[];
	MapIndexed[(i = #2[[1]]; f = #1;
		p = Prod[Pattern[Evaluate[a],BlankNullSequence[]],f[[1]],Pattern[Evaluate[b],BlankNullSequence[]]];
		q = Expand[Evaluate[coeff]*Prod[Evaluate[a],f[[2]],Evaluate[b]]];
		With[{x ={-Evaluate[coeff]*Prod[Evaluate[a]],i,Prod[Evaluate[b]]},y=q},Alternatives[Pattern[Evaluate[coeff],BlankNullSequence[]]*p,p]:> Hold[Sow[x];y]])&
	,sys]//ReleaseHold
]


ExtractRule[f_,i_]:=
Module[{a,b,coeff,p,q},
	a=Unique[];b=Unique[];coeff=Unique[];
	p = Prod[Pattern[Evaluate[a],BlankNullSequence[]],f[[1]],Pattern[Evaluate[b],BlankNullSequence[]]];
	q = Expand[Evaluate[coeff]*Prod[Evaluate[a],f[[2]],Evaluate[b]]];
	With[{x ={-Evaluate[coeff]*Prod[Evaluate[a]],i,Prod[Evaluate[b]]},y=q},Alternatives[Pattern[Evaluate[coeff],BlankNullSequence[]]*p,p]:> Hold[Sow[x];y]]//ReleaseHold
]


SetAttributes[RewriteGroebner,HoldFirst]

RewriteGroebner[cofactors_,F_,OptionsPattern[{Info->False}]]:=
Module[{t,info,N,NC,i,j,k,a,f,b,l,r},
	t = AbsoluteTiming[
	info = OptionValue[Info];
	If[info,
		Print["Rewriting the cofactors has started."];
	];
	N = Length[F];
	NC = Length[cofactors];
	cofactors[[;;N,1,2]] = F;
	
	cofactors = DeleteCases[CollectLeft /@ cofactors,{0,_,_},2];
	
	Monitor[
	Do[
		cofactors[[j]]\[NonBreakingSpace]= Flatten[Reap[
			Do[
				{a,i,b} = cofactors[[j,k]];
				Sow[Map[{Prod[a,#[[1]]],#[[2]],Prod[#[[3]],b]}&,cofactors[[i]]]];
			,{k,Length[cofactors[[j]]]}];
		][[2]],2]//CollectLeft//ExpandLeft;
	,{j,N+1,NC}];
	,NC - j];
	cofactors = ToNonCommutativeMultiply[cofactors];
	][[1]];
	If[info,
		Print["Rewriting the cofactors took in total ", t];
	];
]


(* ::Subsection::Closed:: *)
(*Find Cofactors*)


(* ::Text:: *)
(*ReducedForm[cofactors,G,exp] can be used to reduce the expression exp with the elements of G. The linear combination of these reduction steps is saved in the list cofactors. *)
(*Exp can also be a list of expressions. *)


SetAttributes[ReducedForm,HoldFirst]

ReducedForm[cofactors_,G_,exp_]:=
Module[{t,rules,a,i,b,j},
	cofactors = {};
	rules = ExtractRules[CreateRedSys[G]];
	t = Reap[ToProd[exp]//.rules];
	If[Length[t[[2]]] > 0,
		cofactors = ToNonCommutativeMultiply[Map[({a,j,b} = #; {-a,G[[j]],b})&, t[[2,1]]]];
	];
	ToNonCommutativeMultiply[t[[1]]]
]


ReducedForm[cofactors_,G_,exp:_?ListQ]:=
Module[{i,t,rules,a,b,j},
	cofactors = Table[{},Length[exp]];
	rules = ExtractRules[CreateRedSys[G]];
	Table[
		t = Reap[ToProd[exp[[i]]]//.rules];
		If[Length[t[[2]]] > 0,
			cofactors[[i]] = ToNonCommutativeMultiply[Map[({a,j,b} = #; {-a,G[[j]],b})&, t[[2,1]]]];
		];
		ToNonCommutativeMultiply[t[[1]]]
	,{i,Length[exp]}
	]
]


SetAttributes[ReducedFormIntern,HoldFirst]

ReducedFormIntern[cofactors_,G_,exp:_?ListQ,OptionsPattern[{Shuffle->False}]]:=
Module[{i,t,rules},
	cofactors = Table[{},Length[exp]];
	rules = ExtractRules[CreateRedSys[G]];
	If[OptionValue[Shuffle],
		rules = RandomSample[rules];
	];
	Table[
		t = Reap[ToProd[exp[[i]]]//.rules];
		If[Length[t[[2]]] > 0,
			cofactors[[i]] = ToNonCommutativeMultiply[ReplacePart[#,1->-#[[1]]]&/@t[[2,1]]];
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


GroebnerWithoutCofactors[ideal_,maxiter:_?IntegerQ:10,OptionsPattern[{MaxDeg->Infinity,Info->False,Parallel->True,Sorted->True,Criterion->True}]]:=
Module[{count,spol,p,h,G,lt,info,t,rules,criterion,oldlength,maxdeg,incl,pos,sorted,parallel,syslength,i},

	info = OptionValue[Info];
	criterion = OptionValue[Criterion];
	sorted = OptionValue[Sorted];
	parallel = OptionValue[Parallel];
	maxdeg = OptionValue[MaxDeg];

	G = CreateRedSys[ToProd/@ideal];
	syslength = oldlength = Length[G];
	If[info,Print["G has ", Length[G]," elements in the beginning."];Print[]];
	count = 0; t = 0;

	spol = CheckResolvability2[G,0,Criterion->criterion,MaxDeg->maxdeg,Info->info,Sorted->sorted,Parallel->parallel];
	rules = ExtractRules2[G];

	While[Length[spol] > 0 && count < maxiter,
		i = Length[spol];
		t = AbsoluteTiming[Monitor[Do[
			h = p//.rules;
			If[h =!= 0,  
				lt = LeadingTermIntern[h];
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
	If[OptionValue[Sorted],
		amb = SortBy[amb,Length[#[[1]]]&];
	];
	
	(*generate and reduce S-polynomials*)
	t2 = AbsoluteTiming[
	If[parallel && Length[amb] > 300,
		spol = DeleteDuplicates[DeleteCases[ParallelMap[SPoly2[#,sys[[#[[4,1]]]],sys[[#[[4,2]]]]]//.rules &,amb,DistributedContexts->Automatic,Method->"ItemsPerEvaluation" -> 1000],0]],
		spol = DeleteDuplicates[DeleteCases[Map[SPoly2[#,sys[[#[[4,1]]]],sys[[#[[4,2]]]]]//.rules &,amb],0]]];
	][[1]];

	If[info, Print[Length[spol]," different S-polynomials did not reduce to 0 (computation took ",t2,")"]];
	If[OptionValue[Sorted],
		If[$VersionNumber < 12,
			Sort[spol,SortedQ[LeadingTermIntern[#1][[2]],LeadingTermIntern[#2][[2]]]&],
			SortBy[spol,LeadingTermIntern[#][[2]]&,SortedQ]
		],
		spol
	]
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
	Expand[1/LeadingTermIntern[poly][[1]]*poly]


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
	lt = LeadingTermIntern[poly];
	{lt[[2]], Expand[-1/lt[[1]]*Remainder[poly,lt]]}
]


ToPoly[poly:{_List,_Prod|_Plus|_Times}]:= Prod@@poly[[1]] - poly[[2]]


ToPoly[poly:{_List,0}]:= Prod@@poly[[1]]


ToPoly[sys_]:= ToPoly/@sys;


Rewrite[spolfactors:List[RepeatedNull[List[RepeatedNull[{__,__,__}]]]],cofactors_,OptionsPattern[InputProd->False]]:= 
	Map[Rewrite[#,cofactors]&,spolfactors]


Rewrite[cofactorsFInput_, cofactorsGInput_]:=
Module[{a,b,i,j,l,r,g,rules,occurring,result,G,cofactorsF,cofactorsG},
	If[cofactorsGInput === {},
		cofactorsFInput,
	
		cofactorsF = Map[ToProd,cofactorsFInput,{2}];
		cofactorsG = Map[ToProd,cofactorsGInput,{3}];
		G = ToProd/@(MultiplyOut/@cofactorsGInput);
	
		occurring = Flatten[Position[G,Alternatives@@DeleteDuplicates[cofactorsF[[All,2]]],{1}]];
		rules = Map[
			{a__,G[[#]],b__} -> Sequence@@Map[({l,g,r} = #;{Prod[a,l],g,Prod[r,b]})&,cofactorsG[[#]]]&
		,occurring];
		result = cofactorsF/.rules;

		ToNonCommutativeMultiply[result//CollectLeft//ExpandLeft]
	]
]


SetAttributes[MakeMonic,HoldFirst];

MakeMonic[ideal_]:=
Module[{lc},
	lc = (LeadingTermIntern/@ideal)[[All,1]];
	ideal = Expand[ideal/lc];
	lc
]


Interreduce[ideal_,OptionsPattern[{InputProd->False}]]:=
Module[{G,rules,i,s,gi,cofactors,r,lt,sys,a,b,coeff,p,q,newPart,lc},
	(*set everything up*)
	G = If[OptionValue[InputProd],
		ideal,
		ToProd/@ideal
	];
	lc = MakeMonic[G];
	cofactors = MapIndexed[{{1/#1*Prod[],#2[[1]],Prod[]}}&,lc];
	
	sys = CreateRedSys[ideal];
	a=Unique[];b=Unique[];
	coeff = Unique[];
	rules = ReleaseHold[Table[
		p = Prod[Pattern[Evaluate[a],BlankNullSequence[]],sys[[i,1]],Pattern[Evaluate[b],BlankNullSequence[]]];
		q = Expand[Evaluate[coeff]*Prod[Evaluate[a],sys[[i,2]],Evaluate[b]]];
		With[{x ={-Evaluate[coeff]*Prod[Evaluate[a]],i,Prod[Evaluate[b]]},y=q},Alternatives[Pattern[Evaluate[coeff],BlankNullSequence[]]*p,p]:> Hold[Sow[x];y]]
	,{i,Length[sys]}]];
	i = 1; s = Length[G];
	
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
				lt = LeadingTermIntern[gi];
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


(* ::Subsection:: *)
(*Check Certificate*)


MultiplyOut[certificate:{{_,_,_}...}] := Expand[ToNonCommutativeMultiply[Total[Map[ToProd,certificate]]]]


LinearCombinationQ[certificate:{{_,_,_}...}, assumptions_List] := AllTrue[certificate[[All,2]],MemberQ[assumptions,#]&]


CertificateCoeffQ[_Integer] := True
CertificateCoeffQ[_] := False


CoefficientTest[certificate_] := Module[
	{terms},
	terms = Flatten[Map[ToProd,Flatten[certificate]]/.Plus->List];
	And @@ (MatchQ[0 | Prod[___] | _?CertificateCoeffQ * Prod[___]]/@terms)
];


IntegerCoeffQ[certificate_] := Module[
	{terms},
	terms = Flatten[Map[ToProd,Flatten[certificate]]/.Plus->List];
	And @@ (MatchQ[0 | Prod[___] | _Integer * Prod[___]]/@terms)
];


CheckCertificate[certificate:{{_,_,_}...},claim_, assumptions_] := 
	(MultiplyOut[certificate] === claim) && 
	LinearCombinationQ[certificate,assumptions] && 
	CoefficientTest[certificate]


CheckCertificates[certificates_,claims_,assumptions_] := MapThread[CheckCertificate[#1,#2,assumptions]&,{certificates,claims}]


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


PlotQuiver[Q:Quiver]:= If[$VersionNumber < 12, PlotQuiverOld[Q], PlotQuiverNew[Q]]


PlotQuiverOld[Q_]:= GraphPlot[Map[{#[[2]]->#[[3]],#[[1]]}&,Q],DirectedEdges->True,SelfLoopStyle->0.2]


PlotQuiverNew[Q_]:=
Module[{edgeFun,occured},
	occured = Association@@Map[#[[2;;3]] -> 0 &,Q];
	edgeFun[pts_,e_] := Module[{edge,pos},
		edge = List@@e;
		pos = Position[Q[[All,2;;]],edge,1,++occured[edge]][[-1,1]];
		{Text[Q[[pos,1]],Mean@pts],Arrowheads[0.03],Arrow@pts}
	];
	GraphPlot[
	Map[#[[2]]->#[[3]]&,Q],
		EdgeShapeFunction->edgeFun,DirectedEdges->True,SelfLoopStyle->0.2]
]


TrivialQuiver[F_]:= Module[{vars,v},
	vars = DeleteDuplicates[Select[(ToProd/@F)/.{Prod->List,Plus->List,Times->List}//Flatten,!CoeffQ[#]&]];
	Table[{v,1,1},{v,vars}]
]


QOrderCompatibleQ[p_,Q_]:= Module[{lm},
	lm = ToProd[LeadingTerm[p][[2]]];
	QSignature[lm,Q] === QSignature[p,Q]
]


QConsequenceQ[certificate_,G_,Q_]:=
Module[{f,signaturef,signaturesCertificate,resultsPathCheck,result},
	result = Catch[
	(* check if certificate is indeed a liner combination of elements in F *)
		If[!LinearCombinationQ[certificate,G],
			Print["The input is not a linear combination of elements of the given set."];
			Throw[False]
		];  

		f = MultiplyOut[certificate];
		signaturef = QSignature[f,Q];
	
		(* f has to be compatible *)
		If[signaturef === {}, Throw[False]];
	
		(* Q-consequence test has to pass for (u,v)\in\sigma(f) and every a_i f_i b_i *)
		signaturesCertificate = Map[QSignature[#,Q]&,certificate,{2}];

		(* go through all (u,v)\in\sigma(f) *)
		Do[
			(* given (u,v)\in\sigma(f) check for all a_i f_i b_i *)
			resultsPathCheck = Map[PathCheck[pair,#]&,signaturesCertificate];
			If[MemberQ[resultsPathCheck,False],
				Throw[False]
			];
		,{pair,signaturef}];
		Throw[True];
	];
	result
]


PathCheck[{u_,v_},{siga_,sigg_,sigb_}]:=
Module[{starts,ends,ui,vi},
	starts = Cases[sigb,{u,ui_}->ui];
	ends = Cases[siga,{vi_,v}->vi];
	MemberQ[sigg,{Alternatives@@starts,Alternatives@@ends}]
]


QConsequenceQCrit[certificate_,G_,Q_]:= Module[{f,signaturef},
	(* check if certificate is indeed a liner combination of elements in F *)
	If[!LinearCombinationQ[certificate,G],
		Print["The input is not a linear combination of elements of the given set."];
		Return[False]
	];  

	f = MultiplyOut[certificate];
	signaturef = QSignature[f,Q];
	
	(* f has to be compatible *)
	If[signaturef === {}, Return[False]];
	
	(* Q-consequence criterion has to pass for each summand in the certificate *)
	!MemberQ[Map[QConsequenceCriterion[#,signaturef,Q]&,certificate],False]
]


QConsequenceCriterion[{ai_,gi_,bi_},signaturef_,Q_]:= Module[{sigGi,sigMonomials,mi},
	(* check if there is m_i \in \supp(g_i) s.t. \sigma(m_i) = \sigma(g_i) *)
	sigGi = QSignature[gi,Q];
	sigMonomials = QSignature[MonomialList[gi],Q];
	If[!MemberQ[sigMonomials,sigGi], Return[False]];
	
	(* check if \sigma(f) \subseteq \sigma(a_i m_i b_i) *)
	mi = Extract[MonomialList[gi],FirstPosition[sigMonomials,sigGi,1]];
	SubsetQ[QSignature[Prod[ai,mi,bi],Q],signaturef]
]


(* ::Subsection::Closed:: *)
(*Q-Completion*)


SetAttributes[QCompletion,HoldFirst]

QCompletion[cofactors_,ideal_,Q:Quiver,maxiter:_?IntegerQ:10, OptionsPattern[{Criterion->True,Ignore->0,MaxDeg->Infinity,Info->False,Parallel->True,Sorted->True,IterCount->0}]]:=
Module[{lc,count,spol,lt,info,p,h,G,r,t1,t2,rules,sorted,oldlength,parallel,hrule,maxdeg,intern,criterion,i},
	info = OptionValue[Info];
	sorted = OptionValue[Sorted];
	parallel = OptionValue[Parallel];
	maxdeg = OptionValue[MaxDeg];
	criterion = OptionValue[Criterion];
	intern = FreeQ[ideal,NonCommutativeMultiply]; 

	If[intern,
		G = ideal,
		G = ToProd/@ideal;
		lc = MakeMonic[G];
		cofactors = MapIndexed[{{1/#1*Prod[],#2[[1]],Prod[]}}&,lc];
	];
	
	G = CreateRedSys[G];
	oldlength = Length[G];
	t1 = 0; t2 = 0; count = 0;
	If[info,Print["G has ", Length[G]," elements in the beginning."],Print[]];

	spol = CheckQResolvability[G,Q,OptionValue[Ignore],Criterion->criterion,MaxDeg->maxdeg,Info->info,Sorted->sorted,Parallel->parallel];
	rules = ExtractRules[G];

	While[Length[spol] > 0 && count < maxiter,
		t1 = AbsoluteTiming[
		i = Length[spol];
		Monitor[Do[
			(*reduce it*)
			r = Reap[p[[2]]//.rules]; 
			h = r[[1]];
			(* here we have to check also that h is \leq_Q-compatible and that \sigma(h) \subseteq \sigma(\source(a))*)
			If[h =!= 0 && QOrderCompatibleQ[h,Q] && SubsetQ[QSignature[p[[1]],Q],QSignature[h,Q]],
				If[Length[r[[2]]] > 0,
					p[[3]] = Join[p[[3]],r[[2,1]]]
				];
				lt = LeadingTermIntern[h];
				If[lt[[1]] =!= 1, 
					h = Expand[1/lt[[1]]*h]; p[[3]] = (ReplacePart[#,1 -> 1/lt[[1]]*#[[1]]]&/@ p[[3]])
				];
				hrule = CreateRedSys[h];
				AppendTo[G,hrule];
				AppendTo[cofactors,p[[3]]]; 
				AppendTo[rules,ExtractRule[hrule,Length[G]]];
			];
			i--;
		,{p,spol}];,i];][[1]];
		
		count++;
		If[info, 
			Print["The second reduction took ", t1];
			Print["Iteration ",count + OptionValue[IterCount], " finished. G has now ", Length[G]," elements\n"]
		];
		If[count < maxiter, 
			spol = CheckQResolvability[G,Q,oldlength,Criterion->criterion,MaxDeg->maxdeg,Info->info,Sorted->sorted,Parallel->parallel];
			oldlength = Length[G];
		];
	];
	
	G = ToPoly[G];
	
	If[intern,
		G,
		RewriteGroebner[cofactors,ideal,Info->info];
		ToNonCommutativeMultiply[G]
	]
]


CheckQResolvability[sys_, Q:Quiver, oldlength:_?IntegerQ:0,OptionsPattern[{Criterion->True,MaxDeg->Infinity,Info->False,Parallel->True,Sorted->True}]]:=
Module[{amb,spol,info,rules,parallel,words,r,t1,t2,t3},
	info = OptionValue[Info];
	parallel = OptionValue[Parallel];

	(*generate ambiguities*)
	words = ExtractReducibleWords[sys];
	t1 = AbsoluteTiming[
		amb = GenerateAmbiguities[words[[;;oldlength]],words[[oldlength+1;;]],OptionValue[MaxDeg],Parallel->parallel];
	][[1]];
	
	(* remove ambiguities whose source is not compatible *)
	amb = Select[amb,CompatibleQ[Prod@@#[[1]],Q]&];
	
	If[info,Print[Length[amb]," ambiguities in total (computation took ",t1, ")"]];
	
	(*process ambiguities*)
	If[OptionValue[Criterion],
		amb = DeleteRedundant[amb,words,Info->info];
	];
	If[OptionValue[Sorted],
		amb = SortBy[amb,Length[#[[1]]]&];
	];
	
	(*generate S-polynomials*)
	t2 = AbsoluteTiming[
	spol = DeleteCases[
		If[parallel,
			ParallelMap[QSPoly[#,sys[[#[[4,1]]]],sys[[#[[4,2]]]]]&,amb,DistributedContexts->Automatic,Method->"CoarsestGrained"],
			Map[QSPoly[#,sys[[#[[4,1]]]],sys[[#[[4,2]]]]]&,amb]
		]
	,{_,0,_}];
	][[1]];
	If[info,Print["Generating S-polys: ",t2 ," (",Length[spol], " in total)"]];
	
	(*reduce S-polynomials*)
	rules = ExtractRules[sys];
	(*parallelizing this makes it only slower*)
	t3 = AbsoluteTiming[
	spol = DeleteDuplicatesBy[DeleteCases[Map[(r = Reap[#[[2]]//.rules]; If[r[[1]]=!= 0,
											If[Length[r[[2]]] > 0,
												{#[[1]],r[[1]],Join[#[[3]],r[[2,1]]]},
												{#[[1]],r[[1]],#[[3]]}
												],
											{}
											])&,spol],{}],#[[2]]&];
	][[1]];
	If[info,Print["Reducing S-polys: ",t3, " (",Length[spol], " remaining)"]];
	If[OptionValue[Sorted],
		If[$VersionNumber < 12,
			Sort[spol,SortedQ[LeadingTermIntern[#1[[2]]][[2]],LeadingTermIntern[#2[[2]]][[2]]]&],
			SortBy[spol,LeadingTermIntern[#[[2]]][[2]]&,SortedQ]
		],
		spol
	]
]


QSPoly[amb:_Overlap|_Inclusion,fi_,fj_]:=
Module[{A,C,ABC},
	ABC = Prod@@amb[[1]];
	A = Prod@@amb[[2]];
	C = Prod@@amb[[3]];
	If[amb[[0]]=== Overlap,
			(*Overlap[ABC,A,C]*)
			{ABC, Prod[fi[[2]],C] - Prod[A,fj[[2]]],
				{{A,amb[[4,2]],Prod[]},{-Prod[],amb[[4,1]],C}}},
			(*Inclusion[ABC,A,C)*)
			{ABC, fi[[2]] - Prod[A,fj[[2]],C],
			{{A,amb[[4,2]],C},{-Prod[],amb[[4,1]],Prod[]}}}
	]
]


(* ::Subsection::Closed:: *)
(*Certify*)


Certify[assumptionsInput_List,claimsInput_,Q:Quiver,OptionsPattern[{MaxIter->10,MaxDeg->Infinity,MultiLex->False,Info->False,Parallel->True,Sorted->True,Criterion->True}]]:=
 Module[{info,maxiter,N,shuffle,normalForm,cofactors,cofactorsReduction,G,sigAssump,sigClaim,certificate,toIgnore,toIgnoreOld,alreadyReduced,i,knowns,unknowns,t,assumptions,claims,count},
	info = OptionValue[Info];
	maxiter = OptionValue[MaxIter];
	
	assumptions = ToProd/@assumptionsInput;
	
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

	(*interreduce the generators*)
	{G,cofactors} = Interreduce[assumptions,InputProd->True];
	N = Length[G];
	If[info, Print["\nInterreduced the input from ", Length[assumptions], " polynomials to ", Length[G], ".\n"]];
	
	(*compute the Groebner basis and reduce the claims*)
	If[info, Print["Computing a (partial) Groebner basis and reducing the claim...\n"]];
	(*do computation iteratively*)
	toIgnore = 0;
	i = 1;
	alreadyReduced = Table[False,{i,Length[claims]}];
	certificate = Table[{},{i,Length[claims]}];
	While[MemberQ[alreadyReduced,False] && i <= maxiter,
		toIgnoreOld = Length[G];
		If[info,Print["Starting iteration ", i ,"...\n"]];
		G = Groebner[cofactors,G,1,Ignore->toIgnore,MaxDeg->OptionValue[MaxDeg],Info->OptionValue[Info],Parallel->OptionValue[Parallel],Sorted->OptionValue[Sorted],Criterion->OptionValue[Criterion],IterCount->i-1];
		i++;
		toIgnore = toIgnoreOld;
		count = 0;
		shuffle = False;
		While[MemberQ[alreadyReduced,False] && count < 4,
			count++;
			normalForm = ReducedFormIntern[cofactorsReduction,G,claims,Shuffle->shuffle];
			Do[
				If[alreadyReduced[[j]] === False,
					certificate[[j]] = cofactorsReduction[[j]];
					If[normalForm[[j]] === 0,
						alreadyReduced[[j]] = True;
					];
				];
			,{j,Length[claims]}];
			shuffle = True;
		];
	];
	
	RewriteCertify[certificate,cofactors,assumptionsInput,N];
	certificate = ToNonCommutativeMultiply[Map[ExpandLeft[CollectLeft[#]]&,certificate]];
	
	(* return the reduced claims and the certificates *)
	normalForm = Table[claims[[i]] - MultiplyOut[certificate[[i]]],{i,Length[claims]}];
	If[Head[claimsInput] =!= List, 
		sigClaim = sigClaim[[1]]; certificate = certificate[[1]]; normalForm = normalForm[[1]];
	];
	
	If[info,
		(*full info*)
		If[MemberQ[alreadyReduced,False],
			Print["\nFailed! Not all claims could be reduced to 0."],
			Print["\nDone! All claims were successfully reduced to 0."];
		];
		{Map[QSignature[#,Q]&,assumptions],sigClaim,normalForm,certificate},
		(*only basic info*)
		If[MemberQ[alreadyReduced,False], $Failed, certificate]
	]
]


SetAttributes[RewriteCertify,HoldAll];

RewriteCertify[certificate_,cofactors_,F_,N_]:= 
Module[{a,b,i,j,l,r,f,occurring,done},
	occurring = Select[certificate[[All,All,2]]//Flatten//DeleteDuplicates, # > N&];
	done = {};
	
	(*find all cofactors actually appearing in the certificate*)
	While[Length[occurring] > 0,
		i = occurring[[1]];
		occurring = Delete[occurring,1];
		AppendTo[done,i];
		occurring = Select[Union[occurring,Complement[cofactors[[i,All,2]],done]], # > N&];
	];
	Do[
		cofactors[[i]] = ReplacePart[#,2->F[[#[[2]]]]]&/@cofactors[[i]];
	,{i,N}];
	Do[
		cofactors[[i]] = Map[({a,j,b} = #;Sequence@@Map[({l,f,r} = #;{Prod[a,l],f,Prod[r,b]})&,cofactors[[j]]])&,cofactors[[i]]];
	,{i,Sort[Complement[done,Range[N]]]}];
	Do[
		certificate[[i]] = Map[({a,j,b} = #;Sequence@@Map[({l,f,r} = #;{Prod[a,l],f,Prod[r,b]})&,cofactors[[j]]])&,certificate[[i]]];
	,{i,Length[certificate]}];
]


(* ::Subsection::Closed:: *)
(*End*)


Copyright[a_String,b___String]:= Print[StringJoin[Prepend[{"\n",#}&/@{b},a]]]


Copyright[
    "Package OperatorGB version 1.2.1",
    "Copyright 2019, Institute for Algebra, JKU",
    "by Clemens Hofstadler, clemens.hofstadler@jku.at"];


End[]


EndPackage[]
