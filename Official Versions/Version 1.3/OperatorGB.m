(* ::Package:: *)

(* ::Section::Closed:: *)
(*Documentation*)


(* ::Subsection::Closed:: *)
(*ClearAll*)


BeginPackage["OperatorGB`"]


Global`OperatorGB::usage="";


ClearAll[
	(* Begin *)
	VerboseOperatorGB,
	(* Auxiliary stuff for GB computations *)
	CoeffQ,Prod,ToProd,ToNonCommutativeMultiply,
	Factorize,
	SetUpRing,
	WordOrder,varSets,
	LeadingTerm,DegLex,WeightedDegLex,MultiLex,Weight,SortedQ,
	ExtractReducibleWords,GenerateAmbiguities,Overlap,Inclusion,DeleteRedundant,
	CreateRedSys,ToPoly,Rewrite,
	CollectLeft,CollectRight,ExpandLeft,ExpandRight,
	(* GB computations *)
	Groebner,CheckResolvability,
	GroebnerWithoutCofactors,ApplyRules,
	Intersect,IntersectSubalgebra,ToRightGB,ToLeftGB,IntersectRightIdeal,IntersectLeftIdeal,
	ReducedForm,
	Interreduce,
	Hom,Mon,
	(* Module GB computations *)
	SetModuleBasis,
	ModuleLeadingTerm,
	ModuleAmbiguitiy,ModuleGenerateAmbiguities,
	ModuleGroebner,
	(* Quiver *)
	Quiver,PlotQuiver,
	TrivialQuiver,QuiverFromPolynomials,
	QSignature,CompatibleQ,UniformlyCompatibleQ,QOrderCompatibleQ,QConsequenceQ,QConsequenceQCrit,
	QCompletion,
	(* Certifying operator identities *)
	adj,Pinv,AddAdj,IntegerCoeffQ,
	ApplyLeftCancellability,ApplyRightCancellability,
	Certify,MaxIter,MaxDeg,Parallel,Sorted,Criterion,
	MultiplyOut,LinearCombinationQ,CertificateCoeffQ,IntegerCoeffQ,CheckCertificate,CheckCertificates,
	(* Diagram chasing *)
	ApplyMono,ApplyEpic,ApplyExact,DiagramChase
]


(* ::Subsection::Closed:: *)
(*Begin*)


(*Verbose*)
VerboseOperatorGB::usage="Determines how much information of the computations is printed."


(* ::Subsection::Closed:: *)
(*Auxiliary stuff for GB computations*)


(*Non-commutative multiplication*)
CoeffQ::usage="CoeffQ[k_] should check whether k is in the base ring of the polynomial ring.";
Prod::usage="Prod[m1,m2,...] represents the non-commutative multiplication.";
ToProd::usage="Convert a polynomial from the built in non-commutative multiplication to the Prod data structure";
ToNonCommutativeMultiply::usage="Convert a polynomial from the Prod data structure to the built in non-commutative multiplication";


(*Factorization*)
Factorize::usage="Factorize a single noncommutative polynomial or a list of noncommutative polynomials."


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


(*Ambiguities*)
Overlap::usage="Datastructure for overlap ambiguities."
Inclusion::usage="Datastructure for inclusion ambiguities."
ExtractReducibleWords::usage="Preprocessing before GenerateAmbiguities."
GenerateAmbiguities::usage="GenerateAmbiguities[words] computes all ambigutites among all words in the set 'words'."
DeleteRedundant::usage="Chain criterion to remove redundant ambiguities."


(*Additional stuff*)
CreateRedSys::usage="Converts polynomials into the data structure needed for the Groebner algorithm."
ToPoly::usage="Converts an element of a reduction system back into a polynomial."
Rewrite::usage="Rewrite[cofactorsF, cofactorsG] rewrites a linear combination, which is safed in cofactorsF, with the elements from cofactorsG."
Monomials::usage="Returns a list of all monomials of a given polynomial."


(*Collect cofactors*)
CollectLeft::usage="Tries to collect triples in a cofactor representation having the same left cofactors."
CollectRight::usage="Tries to collect triples in a cofactor representation having the same right cofactors."
ExpandLeft::usage=""
ExpandRight::usage=""


(* ::Subsection::Closed:: *)
(*GB computations*)


(*Groebner basis*)
Groebner::usage="Groebner[cofactors_,ideal_, maxiter:_?IntegerQ:10, OptionsPattern[{Criterion->True,Ignore->0,MaxDeg->Infinity,Parallel->True,Sorted->True}]] executes at most maxiter iterations of the Buchberger algorithm to compute
a (partial) Groebner basis of an ideal. Additionally, for every new element in the Groebner basis a list of cofactors is saved in the list cofactors forming a linear combination of the new element. For further information concerning the OptionPatterns
please see the documentation or the source code."
CheckResolvability::usage=""


(*Groebner basis without cofactors*)
GroebnerWithoutCofactors::usage="GroebnerWithoutCofactors[ideal_,maxiter:_?IntegerQ:10,OptionsPattern[{MaxDeg->Infinity,Parallel->True,Sorted->True,Criterion->True}]] executes at most maxiter iterations 
of the Buchberger algorithm to compute a (partial) Groebner basis of an ideal. For further information concerning the OptionPatterns please see the documentation or the source code."
ApplyRules::usage="ApplyRules[exp,G] reduces the expression exp using polynomials from the set G."


(*Intersection*)
Intersect::usage="Compute a (partial) Groeber basis of the intersection of two two-sided ideals."
IntersectSubalgebra::usage="Compute a (partial) Groebner basis of the intersection of a two-sided ideal with a (finitely generated) subalgebra."
ToRightGB::usage="Converts a two-sided Gr\[ODoubleDot]bner basis of an ideal into a right Gr\[ODoubleDot]bner basis."
ToLeftGB::usage="Converts a two-sided Gr\[ODoubleDot]bner basis of an ideal into a left Gr\[ODoubleDot]bner basis."
IntersectRightIdeal::usage="Computes a (partial) right Gr\[ODoubleDot]bner basis of the intersection of a two-sided ideal with a right sided ideal."
IntersectLeftIdeal::usage="Computes a (partial) left Gr\[ODoubleDot]bner basis of the intersection of a two-sided ideal with a left sided ideal."


(*Computing normal forms*)
ReducedForm::usage="ReducedForm[cofactors,G,exp] reduces the expression exp by the elements of G and saves the cofactors of the reduction process in the list cofactors. The argument exp can also
be a list of expressions, then all expressions are reduced."


(*Interreduction*)
Interreduce::usage="Interreduces a set of polynomials."


(*Homogeneous & Monomial part*)
Hom::usage="Computes a generating set of the ideal generated by all homogeneous polynomials in a two-sided ideal."
Mon::usage="Computes a generating set of the ideal generated by all monomial in a right ideal."


(* ::Subsection::Closed:: *)
(*Module GB computations*)


(* Module basis *)
SetModuleBasis::usage="Define the basis elements of the module."


(* Module Ordering *)
ModuleLeadingTerm::usage="Leading term of the polynomial w.r.t. the specified module ordering."


(* Module Ambiguities *)
ModuleAmbiguitiy::usage="Data structure of module S-polynomials."
ModuleGenerateAmbiguities::usage="Computes all module ambigutites among all words in the set 'words'."


(* Module Groebner basis *)
ModuleGroebner::usage="ModuleGroebner[cofactors_,M_, maxiter:_?IntegerQ:10, OptionsPattern[{Ignore->0,MaxDeg->Infinity,Parallel->True,Sorted->True}]] executes at most maxiter iterations of the Buchberger algorithm to compute
a (partial) Groebner basis of a K<X>-bimodule. Additionally, for every new element in the Groebner basis a list of cofactors is saved in the list cofactors forming a linear combination of the new element. For further information concerning the OptionPatterns
please see the documentation or the source code."


(* ::Subsection::Closed:: *)
(*Quiver*)


(*Quiver data structure & basic functionality*)
Quiver::usage="Data structure of a Quiver"
PlotQuiver::usage="Plot a quiver Q"


(*Constructing quivers from polynomials*)
TrivialQuiver::usage="Returns the trivial quiver containing all variables of the given input."
QuiverFromPolynomials::usage="Returns the most general quiver (with unique labels) such a given set of polynomials is (uniformly) compatible."


(*Compatibility checks*)
QSignature::usage="QSignature[poly,Q] returns the signature of the polynomial poly w.r.t. the quiver Q (not necessarily with unique labels)"
CompatibleQ::usage="Tests whether a polynomial is compatible with a quiver."
UniformlyCompatibleQ::usage="Tests whether a polynomial is uniformly compatible with a quiver."
QOrderCompatibleQ::usage="Tests whether a polynomial is Q-order-compatible with a quiver and the order defined by SetUpRing."
QConsequenceQ::usage="Tests whether a certificate is a Q-consequence of a set of polynomials and a quiver using the definition of Q-consequence."
QConsequenceQCrit::usage="Tests whether a certificate is a Q-consequence of a set of polynomials and a quiver using a criterion."


(*Q-completion*)
QCompletion::usage="Q-completion procedure"


(* ::Subsection::Closed:: *)
(*Certifying operator identities*)


(*Adjungate operator definition & auxiliary stuff*)
adj::usage="adj[A] represents the adjoint of the operator A"
Pinv::usage="Gives the 4 Moore-Penrose equations."
AddAdj::usage="Adds the adjoint statements to a given list of statements."
IntegerCoeffQ::usage="Tests if a certficate contains only integer coefficients."


(*Cancellability*)
ApplyLeftCancellability::usage="Tries to find elements in a given Ideal I to which left cancellability of a given element can applied."
ApplyRightCancellability::usage="Tries to find elements in a given Ideal I to which right cancellability of a given element can applied."


(*Certify*)
Certify::usage="Certifies whether a certain claim is a consequence of some assumptions via Groebner 
basis computations. Additionally, compatibility with a given quiver is checked."


(*Check certificate*)
MultiplyOut::usage="To multiply out a list of cofactors given in terms of the built in non-commutative multiplication."
LinearCombinationQ::usage="Checks whether a given set of triples is a linear combination of a set of polynomials."
CertificateCoeffQ::usage=""
IntegerCoeffQ::usage="Checks whether a given certificate only contains integer coefficients"
CheckCertificate::usage="Check that a single certificate indeed gives the claim, is w.r.t. to the assumptions and that it only consists of integer coefficients"
CheckCertificates::usage="Applies CheckCertificate to several certificates"


(* ::Subsection::Closed:: *)
(*Diagram chasing*)


ApplyMono::usage=""
ApplyEpic::usage=""
ApplyExact::usage=""
DiagramChase::usage="Does an automated diagram chase."


(* ::Section::Closed:: *)
(*Begin*)


Begin["`Private`"]


VerboseOperatorGB = 0


(* ::Section::Closed:: *)
(*Auxiliary stuff for GB computations*)


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


ToNonCommutativeMultiply[poly_]:= poly//.{Prod[]:>1, Prod[a_]:>Times[a], Prod[a_,b__]:>NonCommutativeMultiply[a,b]}


(* ::Subsection::Closed:: *)
(*Factorization*)


FProd[a___,FProd[b___],c___]:=FProd[a,b,c]
FProd[FProd[a___]-FProd[b___]]:=FProd[a]-FProd[b]
FProd[FProd[a___]+FProd[b___]]:=FProd[a]+FProd[b]


Factorize[FF_List]:=Module[{F,l1,l2,l3,r1,r2,r3,a,b,l,r,x,y,c,c1,c2},
	F = (ToProd/@FF)/.Prod->FProd;
	
	l1 = Plus[l___,FProd[x__,a___],FProd[x__,b___],r___]:>Plus[l,FProd[x,FProd[a]+FProd[b]],r];
	l2 = Plus[l___,FProd[x__,a___],Times[c_,FProd[x__,b___]],r___]:>Plus[l,FProd[x,FProd[a]+c*FProd[b]],r];
	l3 = Plus[l___,Times[c1_,FProd[x__,a___]],Times[c2_,FProd[x__,b___]],r___]:>Plus[l,FProd[x,c1*FProd[a]+c2*FProd[b]],r];

	r1 = Plus[l___,FProd[a___,y__],FProd[b___,y__],r___]:>Plus[l,FProd[FProd[a]+FProd[b],y],r];
	r2 = Plus[l___,FProd[a___,y__],Times[c_,FProd[b___,y__]],r___]:>Plus[l,FProd[FProd[a]+c*FProd[b],y],r];
	r3 = Plus[l___,Times[c1_,FProd[a___,y__]],Times[c2_,FProd[b___,y__]],r___]:>Plus[l,FProd[c1*FProd[a]+c2*FProd[b],y],r];
	
	(F//.{l1,l2,l3,r1,r2,r3})//.{FProd[]:>1, FProd[a_]:>Times[a], FProd[a_,b__]:>NonCommutativeMultiply[a,b]}
]


Factorize[ff_]:=Module[{f,l1,l2,l3,r1,r2,r3,a,b,l,r,x,y,c,c1,c2},
	f = ToProd[ff]/.Prod->FProd;
	
	l1 = Plus[l___,FProd[x__,a___],FProd[x__,b___],r___]:>Plus[l,FProd[x,FProd[a]+FProd[b]],r];
	l2 = Plus[l___,FProd[x__,a___],Times[c_,FProd[x__,b___]],r___]:>Plus[l,FProd[x,FProd[a]+c*FProd[b]],r];
	l3 = Plus[l___,Times[c1_,FProd[x__,a___]],Times[c2_,FProd[x__,b___]],r___]:>Plus[l,FProd[x,c1*FProd[a]+c2*FProd[b]],r];

	r1 = Plus[l___,FProd[a___,y__],FProd[b___,y__],r___]:>Plus[l,FProd[FProd[a]+FProd[b],y],r];
	r2 = Plus[l___,FProd[a___,y__],Times[c_,FProd[b___,y__]],r___]:>Plus[l,FProd[FProd[a]+c*FProd[b],y],r];
	r3 = Plus[l___,Times[c1_,FProd[a___,y__]],Times[c2_,FProd[b___,y__]],r___]:>Plus[l,FProd[c1*FProd[a]+c2*FProd[b],y],r];
	
	(f//.{l1,l2,l3,r1,r2,r3})//.{FProd[]:>1, FProd[a_]:>Times[a], FProd[a_,b__]:>NonCommutativeMultiply[a,b]}
]


(* ::Subsection::Closed:: *)
(*Setting up the ring*)


SetUpRing[vars_List]:= (
	WordOrder = vars;
	SortedQ := DegLex;
	If[VerboseOperatorGB > 1,
		Print[Sequence@@Map[ToString[#,StandardForm]<>" < "&,vars[[;;-2]]] <> ToString[vars[[-1]],StandardForm]];
	];
);


SetUpRing[S1_List,SRest:_List..]:= Module[{string},
	WordOrder = Flatten[Join[S1,SRest]]; 
	varSets = {S1,SRest};
	SortedQ := MultiLex;
	If[VerboseOperatorGB > 1,
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


DeleteRedundantSimple[amb_List,lt_List]:= 
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
	If[VerboseOperatorGB > 1,
		Print["Removed ", Length[amb] - Length[result], " ambiguities in ",t]];
	result
]


DeleteRedundantComplex[amb_List,lt_List]:= 
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
	If[VerboseOperatorGB > 1,
		Print["Removed ", Length[amb] - Length[result], " ambiguities in ",t]];
	result
]


DeleteRedundant[amb_List,lt_List]:= 
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
	If[VerboseOperatorGB > 1,
		Print["Removed ", Length[amb] - Length[result], " ambiguities in ",t]];
	result
]


(* ::Text:: *)
(*Gebauer-M\[ODoubleDot]ller criterion*)


GebauerMoeller[ambInput_List,lt_List]:= 
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
	If[VerboseOperatorGB > 1,
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


Monomials[f_]:=Module[{c,x},
	MonomialList[f]/.{c__*Prod[x___]:>{x},Prod->List}
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


(* ::Section::Closed:: *)
(*GB computations*)


(* ::Subsection::Closed:: *)
(*Groebner basis*)


(* ::Text:: *)
(*Implementation of the Buchberger algorithm to compute a (partial) Groebner basis of the ideal 'ideal' with at most 'maxiter' iterations being executed (default: 10). The return value is a  set of polynomials G = {f1,...,fn,g1,...gm} consisting of the elements f1,....,fn from 'ideal' and new elements g1,...,gm. For every element g\in G a list forming a linear combination of g consisting of elements from 'ideal' and certain cofactors is saved in the list cofactors.*)
(**)
(*OptionPattern:*)
(*	- Criterion (default: True): Tries to detect and delete redundant ambiguities during the Groebner basis computation.*)
(*	- Ignore (default: 0): A non-negative integer that determines how many elements of the input will be ignored during the first computation of the ambiguities. *)
(*	- MaxDeg (default: Infinity): Only ambiguities with degree smaller than or equal to MaxDeg will be considered during the Groebner basis computation (larger ambiguities are simply ignored). *)
(*	- Parallel (default: True): Determines whether the computations for which it is possible, are executed in parallel (which speeds up the computation) or in series.*)
(*	- Sorted (default: True):  Sorts the ambiguities before processing in ascending order. This speeds up the computation but results in a different (partial) Groebner basis.*)
(*	- IterCount (default: 0): defines from which number the iterations are counted (only relevant for the printed information)*)
(**)


SetAttributes[Groebner,HoldFirst]

Groebner[cofactors_,ideal_, maxiter:_?IntegerQ:10, OptionsPattern[{Criterion->True,Ignore->0,MaxDeg->Infinity,Parallel->True,Sorted->True,IterCount->0}]]:=
Module[{lc,count,spol,lt,p,h,G,r,t1,t2,rules,sorted,oldlength,parallel,hrule,maxdeg,intern,criterion,i,itercount},
	sorted = OptionValue[Sorted];
	parallel = OptionValue[Parallel];
	maxdeg = OptionValue[MaxDeg];
	criterion = OptionValue[Criterion];
	itercount = OptionValue[IterCount];
	intern = !FreeQ[ideal,Prod]; 
	
	If[Head[cofactors] =!= List,
		cofactors = {};
	];

	If[intern,
		G = ideal,
		G = ToProd/@ideal;
		lc = MakeMonic[G];
		cofactors = MapIndexed[{{1/#1*Prod[],#2[[1]],Prod[]}}&,lc];
	];
	
	G = CreateRedSys[G];
	oldlength = OptionValue[Ignore];
	t1 = 0; t2 = 0; count = 0;
	If[VerboseOperatorGB > 0,
		Print["G has ", Length[G]," elements in the beginning.\n"]
	];

	rules = ExtractRules[G];
	spol = {0};

	While[Length[spol] > 0 && count < maxiter,
		count++;
		If[VerboseOperatorGB > 0,
			Print["Starting iteration ", count + itercount ,"..."]
		];
		spol = CheckResolvability[G,oldlength,Criterion->criterion,MaxDeg->maxdeg,Sorted->sorted,Parallel->parallel];
		oldlength = Length[G];
		
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
		
		If[VerboseOperatorGB > 1, 
			Print["The second reduction took ", t1]
		];
		If[VerboseOperatorGB > 0,
			Print["Iteration ",count + itercount, " finished. G has now ", Length[G]," elements\n"]
		];
	];
	
	G = ToPoly[G];
	
	If[intern,
		G,
		RewriteGroebner[cofactors,ideal];
		ToNonCommutativeMultiply[G]
	]
]


(* ::Text:: *)
(*CheckResolvability returns all S-polynomials from the reduction system sys which can not be reduced to zero. Additionally, for each S-polynomial a list containing the linear combination how the S-polynomial was generated from the elements of sys is returned.*)
(*For a description of the OptionPatterns see the documentation of the Groebner method.*)


CheckResolvability[sys_,oldlength:_?IntegerQ:0,OptionsPattern[{Criterion->True,MaxDeg->Infinity,Parallel->True,Sorted->True}]]:=
Module[{amb,spol,rules,parallel,words,r,t1,t2,t3,str},
	parallel = OptionValue[Parallel];

	(*generate ambiguities*)
	words = ExtractReducibleWords[sys];
	t1 = AbsoluteTiming[
		amb = GenerateAmbiguities[words[[;;oldlength]],words[[oldlength+1;;]],OptionValue[MaxDeg],Parallel->parallel];
	][[1]];
	If[VerboseOperatorGB > 0,
		str = ToString[Length[amb]] <> " ambiguities in total" <> If[VerboseOperatorGB > 1, " (computation took " <> ToString[t1] <> ")",""];
		Print[str];
	];
	
	(*process ambiguities*)
	If[OptionValue[Criterion],
		amb = DeleteRedundant[amb,words];
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
	
	If[VerboseOperatorGB > 1,
		Print["Generating S-polys: ",t2 ," (",Length[spol], " in total)"]
	];
	
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
	If[VerboseOperatorGB > 1,
		Print["Reducing S-polys: ",t3, " (",Length[spol], " remaining)"]
	];
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

RewriteGroebner[cofactors_,F_]:=
Module[{t,N,NC,i,j,k,a,b},
	t = AbsoluteTiming[
	If[VerboseOperatorGB > 0,
		Print["Rewriting the cofactors..."];
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
	If[VerboseOperatorGB > 1,
		Print["Rewriting the cofactors took in total ", t];
	];
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
(*	- Parallel (default: True): Determines whether the computations for which it is possible, are executed in parallel (which speeds up the computation) or in series.*)
(*	- Sorted (default: True):  Sorts the ambiguities before processing in ascending order. This speeds up the computation but results in a different (partial) Groebner basis.*)
(*	*)


GroebnerWithoutCofactors[ideal_,maxiter:_?IntegerQ:10,OptionsPattern[{MaxDeg->Infinity,Parallel->True,Sorted->True,Criterion->True}]]:=
Module[{count,spol,p,h,G,lt,t,rules,criterion,oldlength,maxdeg,sorted,parallel,i},

	criterion = OptionValue[Criterion];
	sorted = OptionValue[Sorted];
	parallel = OptionValue[Parallel];
	maxdeg = OptionValue[MaxDeg];

	G = CreateRedSys[ToProd/@ideal];
	If[VerboseOperatorGB > 0,
		Print["G has ", Length[G]," elements in the beginning.\n"]
	];
	
	count = 0; t = 0; oldlength = 0;
	rules = ExtractRules2[G];
	spol = {0};

	While[Length[spol] > 0 && count < maxiter,
		count++;
		If[VerboseOperatorGB > 0,
			Print["Starting iteration ", count,"..."]
		];
		spol = CheckResolvability2[G,oldlength,Criterion->criterion,MaxDeg->maxdeg,Sorted->sorted,Parallel->parallel];
		oldlength = Length[G];
		
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
		If[VerboseOperatorGB > 1, 
			Print["The second reduction took ", t]
		];
		If[VerboseOperatorGB > 0,
			Print["Iteration ",count, " finished. G has now ", Length[G]," elements\n"]
		];
	];
	
	ToNonCommutativeMultiply[ToPoly[G]]
]


CheckResolvability2[sys_,oldlength:_?IntegerQ:0,OptionsPattern[{Criterion->True,MaxDeg->Infinity,Sorted->True,Parallel->True}]]:=
Module[{amb,spol,t1,t2,lists,rules,words,parallel,str},
	parallel = OptionValue[Parallel];

	words = ExtractReducibleWords[sys];
	rules = ExtractRules2[sys];
	
	(*generate ambiguities*)
	t1 = AbsoluteTiming[
		amb = GenerateAmbiguities[words[[;;oldlength]],words[[oldlength+1;;]],OptionValue[MaxDeg],Parallel->True]
	][[1]];
	If[VerboseOperatorGB > 0,
		str = ToString[Length[amb]] <> " ambiguities in total" <> If[VerboseOperatorGB > 1, " (computation took " <> ToString[t1] <> ")", ""];
		Print[str];
	];
	
	(*process ambiguities*)
	If[OptionValue[Criterion],
		amb = DeleteRedundant[amb,words[[All,1]]]
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

	If[VerboseOperatorGB > 1,
		Print["Reducing S-polys: ", t2, " (", Length[spol], "remaining)"]
	];
	
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
(*Intersections*)


(* ::Subsubsection::Closed:: *)
(*with two-sided ideal*)


Intersect[II_,JJ_,OptionsPattern[{Vars->{},MaxIter->1, CutOff->-1}]]:=
Module[{t,I,J,IJ,oldIsDegLex,cofactors,G,vars,spol,rules,commutator,a,b,c,t1,t2,N,pos,F,cutoff,red},
	cutoff = OptionValue[CutOff];
	
	(* convert to Prod data structure *)
	If[FreeQ[II,Prod]||FreeQ[JJ,Prod],
		I = ToProd/@II; J = ToProd/@JJ,
		I = II; J = JJ
	];
	
	(* consider ideals in K<vars>  *)
	vars = If[Length[OptionValue[Vars]] > 0,
			OptionValue[Vars],
			WordOrder
			];
	
	(* set up elimination order *)
	If[VerboseOperatorGB > 1,
		Print["Switching to elimination order..."];
	];
	oldIsDegLex = SortedQ === DegLex;
	If[oldIsDegLex,
		SetUpRing[WordOrder,{t}],
		SetUpRing[Sequence@@Reverse[varSets],{t}]
	];
	
	(* set up the ideal (tI + (1-t)J + tX - Xt) *)
	G = Join[Map[Prod[t,#] &,I],Map[Prod[1-t,#] &,J]];
	commutator = Prod[a___,t,b_,c___]->Prod[a,b,t,c];
	
	(* compute Groebner basis *)
	N = 0;
	Do[
		If[VerboseOperatorGB > 0,
			Print["\nStarting iteration ", i, "..."]
		];
		(* do interreduction *)
		(* but only use rules until cutoff - speeds things up *)
		If[0 < cutoff && cutoff < Length[G],
			F = Join[G[[;;cutoff]],Select[G[[cutoff+1;;]],FreeQ[#,t]&]];
			rules = Append[ExtractRules[CreateRedSys[F]],commutator];
			If[VerboseOperatorGB > 0,
				Print["No interreduction due to cutoff."];
			],
			t1 = AbsoluteTiming[
				rules = Append[ExtractRules[CreateRedSys[G]],commutator];
				{G,cofactors} = Interreduce[G,Rules->rules];
				rules = Append[ExtractRules[CreateRedSys[G]],commutator];
			][[1]];
			If[VerboseOperatorGB > 1,
				Print["Interreduction took ", t1, ". G has now ", Length[G], " elements."];
			];
		];
		
		(* generate S-polies *)
		(* only pick those elements in G that have not already been there before *)
		pos = Position[cofactors,{{_,i_,_}}/;i <= N];
		pos = Complement[List/@Range[Length[G]],pos];
		spol = SPolyForIntersection[Extract[G,pos],t,vars];
		
		(* reduce S-polies - do this in blocks of smaller size - speeds things up  *)
		t2 = AbsoluteTiming[
			red = {};
			Do[
				red = Join[red,DeleteCases[spol[[1000*(j-1) + 1 ;; Min[1000*j,Length[spol]]]]//.rules,0]]//DeleteDuplicates;
			,{j,Length[spol]/1000+1}]
		][[1]];
		
		N = Length[G];
		G = Join[G,red];
		If[VerboseOperatorGB > 1,
			Print["Reduction took ", t2 , " (", Length[red], " reminaing)"];
			Print["G has now ", Length[G], " elements (", Count[G,x_/;FreeQ[x,t]&], " not involving t)"];
		];
	,{i,OptionValue[MaxIter]}];
	
	(* set order back to what it was before *)
	If[VerboseOperatorGB > 1,
		Print["\nSwitching back to original order..."];
	];
	If[oldIsDegLex,
		SetUpRing[WordOrder[[;;-2]]],
		SetUpRing[Sequence@@Reverse[varSets[[2;;]]]]
	];
	
	(* return all elements not involving t *)
	ToNonCommutativeMultiply/@Select[G,FreeQ[#,t]&]
]


SPolyForIntersection[G_,t_,vars_]:=
Module[{mon,f,g,lm,tail,spol},
	(* generate new S-polies *)
	mon = Map[MonomialList,G];
	(* every element is of the form ft + g *)
	f = (Select[#,!FreeQ[#,t]&]&/@ mon)/.t->1;
	f = Apply[Plus,f,{1}];
	g = Select[#,FreeQ[#,t]&]&/@ mon;
	g = Apply[Plus,g,{1}];

	lm = Prod@@LeadingTermIntern[#][[2]]&/@f;
	tail = f-lm;
	spol = (helper[lm,tail,g,t,#]&/@ vars)//Flatten;
	If[VerboseOperatorGB > 0,
		Print[Length[spol], " S-polies"];
	];
	spol
]

helper[lm_,tail_,g_,t_,x_] := MapThread[Prod[#1,x,t] + Prod[#2,t,x] + Prod[#3,x]&,{lm,tail,g}]


(* ::Subsubsection::Closed:: *)
(*with subalgebra*)


IntersectSubalgebra[I_,subalg_,OptionsPattern[{MaxIter->5, SplitAt->-1}]]:=
Module[{G,Y,IJ,cofactors,oldIsDegLex,splitAt},
	splitAt = OptionValue[SplitAt];

	
	(* set up generating set {yi - gi, fj |\[NonBreakingSpace]gi \in subalg, fj \in I} *)
	Y = Map[Unique[]&,subalg];
	
	IJ = ToProd/@Join[Y[[;;splitAt]] - subalg[[;;splitAt]], I, Y[[splitAt+1;;]] - subalg[[splitAt+1;;]]];
	
	(* switch to elimination order *)
	oldIsDegLex = SortedQ === DegLex;
	If[VerboseOperatorGB > 1,
		Print["Switching to elimination order..."];
	];
	If[oldIsDegLex,
		SetUpRing[Y,WordOrder],
		SetUpRing[Y,Sequence@@Reverse[varSets]]
	];
	
	(* compute Groebner basis *)
	cofactors = {};
	G = Groebner[cofactors,IJ,OptionValue[MaxIter]];
	If[VerboseOperatorGB > 1,
		Print["\nSwitching back to original order..."];
	];
	
	(* switch back to old order *)
	If[oldIsDegLex,
		SetUpRing[WordOrder[[Length[Y]+1;;]]],
		SetUpRing[Sequence@@Reverse[varSets[[;;-2]]]]
	];
	
	(* 1 \in G => I \cap subalg = subalg *)
	If[MemberQ[G,Prod[]],
		ToNonCommutativeMultiply/@subalg,
		(* pick all non-zero elements in K<Y> and replace the Subscript[y, i] *)
		ToNonCommutativeMultiply/@DeleteCases[Select[G,FreeQ[#,Alternatives@@WordOrder]&]/.Thread[Rule[Y,ToProd/@subalg]],0]
	]
]


(* ::Subsubsection::Closed:: *)
(*with right ideal*)


ToRightGB[GG_,Q_,vvars_,length_]:=
Module[{G,targetG,G2,lmAll,lmPrefix,m,vars,new,sourcem,t0,t1,t2,intern},
	intern = !FreeQ[GG,Prod];
	
	(* get a basis for (lm(G)) *)
	G = Interreduce[GG][[1]];
	lmAll = (LeadingTerm/@G)[[All,2]];
	targetG = (QSignature[#,Q]&/@Prod/@lmAll)[[All,All,2]];
	lmPrefix = Alternatives@@ Map[Prod[___,Sequence@@#,___]&,lmAll];
	lmAll = Alternatives@@ Map[Prod[___,Sequence@@#,__]&,lmAll];
	G = Thread[List[G,targetG]];
	
	(* pick those variables that are not in (lm(G)) *)
	vars = DeleteCases[ToProd/@vvars,lmPrefix];
	m = {{1,1}};
	
	G2 = Reap[Do[
		m = DeleteCases[Outer[Prod,vars,m[[All,1]]]//Flatten,lmPrefix];
		sourcem = (QSignature[#,Q]&/@m)[[All,All,1]];
		m = Thread[List[m,sourcem]];
		m = DeleteCases[m,{_,{}}];
		new = Reap[
			Outer[
				If[ContainsAny[#1[[2]],#2[[2]]],
					Sow[Prod[#1[[1]],#2[[1]]]],
					Sequence[]
				]&
			,m,G,1,1];
		][[2]]//Flatten;
		new = DeleteCases[new,lmAll];
		If[Length[new] === 0,
			If[VerboseOperatorGB > 0,
				Print["All compatible polynomials found at iteration ",i];
			];
			Break[],
		Sow[new]
		];
	,{i,length}];
	][[2]]//Flatten;
	G = Join[G[[All,1]],G2];
	If[intern,
		G,
		ToNonCommutativeMultiply/@G
	]
]


IntersectRightIdeal[GG_,RR_,Q_,OptionsPattern[{Vars->{},Length->5}]]:=
Module[{vars,G,R,rightGB,t,int,oldIsDegLex},
	(* fix variables *)
	vars = OptionValue[Vars];
	If[Length[vars]===0,
		vars = WordOrder
	];

	(* convert to Prod data structure *)
	G = ToProd/@GG;
	R = ToProd/@RR;

	(* compute a right compatible GB of G *)
	rightGB = ToRightGB[G,Q,vars,OptionValue[Length]];

	(* set up elimination order *)
	t = Unique[];
	If[VerboseOperatorGB > 1,
		Print["Switching to elimination order..."];
	];
	oldIsDegLex = SortedQ === DegLex;
	If[oldIsDegLex,
		SetUpRing[WordOrder,{t}],
		SetUpRing[Sequence@@Reverse[varSets],{t}]
	];

	(* intersect the two right ideals *)
	int = Join[Prod[t,#]&/@rightGB,Prod[1-t,#]&/@R];
	int = Interreduce[int,OneSided->"right"][[1]];
	
	(* set order back to what it was before *)
	If[VerboseOperatorGB > 1,
		Print["\nSwitching back to original order..."];
	];
	If[oldIsDegLex,
		SetUpRing[WordOrder[[;;-2]]],
		SetUpRing[Sequence@@Reverse[varSets[[2;;]]]]
	];

	ToNonCommutativeMultiply/@Select[int,FreeQ[#,t]&]
]


(* ::Subsubsection::Closed:: *)
(*with left ideal*)


ToLeftGB[GG_,Q_,vvars_,length_]:=
Module[{G,G2,lmAll,lmPostfix,m,vars,new,intern},
	intern = !FreeQ[GG,Prod];
	(* get a basis for (lm(G)) *)
	G = Interreduce[GG][[1]];
	lmAll = (LeadingTerm/@G)[[All,2]];
	lmPostfix = Alternatives@@ Map[Prod[___,Sequence@@#,___]&,lmAll];
	lmAll = Alternatives@@ Map[Prod[__,Sequence@@#,___]&,lmAll];

	(* pick those variables that are not in (lm(G)) *)
	vars = DeleteCases[ToProd/@vvars,lmPostfix];
	m = {1};
	G2 = Reap[Do[
		m = DeleteCases[Outer[Prod,m,vars]//Flatten,lmPostfix];
		m = Select[m,CompatibleQ[#,Q]&];
		new = DeleteCases[Outer[Prod,G,m]//Flatten,lmAll];
		If[Length[new] === 0,
			If[VerboseOperatorGB > 0,
				Print["All compatible polynomials found at iteration ",i];
			];
			Break[],
		Sow[new]
		];
	,{i,length}];
	][[2]]//Flatten;
	G = Join[G,Select[G2,CompatibleQ[#,Q]&]];
	If[intern,
		G,
		ToNonCommutativeMultiply/@G
	]
]


IntersectLeftIdeal[GG_,LL_,Q_,OptionsPattern[{Vars->{},Length->5}]]:=
Module[{vars,G,L,leftGB,t,int,oldIsDegLex},
	(* fix variables *)
	vars = OptionValue[Vars];
	If[Length[vars]===0,
		vars = WordOrder
	];

	(* convert to Prod data structure *)
	G = ToProd/@GG;
	L = ToProd/@LL;

	(* compute a left compatible GB of G *)
	leftGB = ToLeftGB[G,Q,vars,OptionValue[Length]];

	(* set up elimination order *)
	t = Unique[];
	If[VerboseOperatorGB > 1,
		Print["Switching to elimination order..."];
	];
	oldIsDegLex = SortedQ === DegLex;
	If[oldIsDegLex,
		SetUpRing[WordOrder,{t}],
		SetUpRing[Sequence@@Reverse[varSets],{t}]
	];

	(* intersect the two left ideals *)
	int = Join[Prod[#,t]&/@leftGB,Prod[#,1-t]&/@L];
	int = Interreduce[int,OneSided->"left"][[1]];

	(* set order back to what it was before *)
	If[VerboseOperatorGB > 1,
		Print["\nSwitching back to original order..."];
	];
	If[oldIsDegLex,
		SetUpRing[WordOrder[[;;-2]]],
		SetUpRing[Sequence@@Reverse[varSets[[2;;]]]]
	];

	ToNonCommutativeMultiply/@Select[int,FreeQ[#,t]&]
]


(* ::Subsection::Closed:: *)
(*Computing normal forms*)


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
(*Interreduction*)


Interreduce[ideal_,OptionsPattern[{Rules->None,Complete->True,OneSided->False}]]:=
If[OptionValue[Complete],
	interreduceFull[ideal,Rules->OptionValue[Rules],OneSided->OptionValue[OneSided]],
	interreduceHead[ideal,Rules->OptionValue[Rules],OneSided->OptionValue[OneSided]]
]


interreduceFull[ideal_,OptionsPattern[{Rules->None,OneSided->False}]]:=
Module[{G,rules,i,s,gi,cofactors,r,lt,sys,a,b,coeff,p,q,newPart,lc,oneSided},
	oneSided = ToLowerCase[OptionValue[OneSided]];
	(*set everything up*)
	G = If[FreeQ[ideal,Prod],
		ToProd/@ideal,
		ideal
	];
	lc = MakeMonic[G];
	cofactors = MapIndexed[{{1/#1*Prod[],#2[[1]],Prod[]}}&,lc];
	
	sys = CreateRedSys[ideal];
	rules = If[OptionValue[Rules] === None,
		Which[oneSided === "left",
				ExtractLeftRules[sys],
			oneSided === "right",
				ExtractRightRules[sys],
			True,
				ExtractRules[sys]
		],
		OptionValue[Rules]
	];
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
				(* if rules was used for tracing, adapt cofactors *)
				If[Length[r[[2]]] > 0,
					r = r[[2,1]];
					newPart = Flatten[Table[Map[{Prod[r[[i,1]],#[[1]]],#[[2]],Prod[#[[3]],r[[i,3]]]}&,cofactors[[r[[i,2]]]]],{i,Length[r]}],1];
					cofactors[[i]] = Join[cofactors[[i]],newPart];
				];
				(* make poly monic if it is not already *)
				If[lt[[1]] =!= 1, 
					gi = Expand[1/lt[[1]]*gi]; 
					cofactors[[i]] = (ReplacePart[#,1 -> 1/lt[[1]]*#[[1]]]&/@ cofactors[[i]]);
				];
				G[[i]]=gi;
				sys = CreateRedSys[gi];
				Which[oneSided === "left",
						p = Prod[Pattern[Evaluate[a],BlankNullSequence[]],sys[[1]]];
						q = Expand[Evaluate[coeff]*Prod[Evaluate[a],sys[[2]]]];
						rules[[i]]= ReleaseHold[With[{x ={-Evaluate[coeff]*Prod[Evaluate[a]],i,Prod[]},y=q},Alternatives[Pattern[Evaluate[coeff],BlankNullSequence[]]*p,p]:> Hold[Sow[x];y]]],
					oneSided === "right",
						p = Prod[sys[[1]],Pattern[Evaluate[b],BlankNullSequence[]]];
						q = Expand[Evaluate[coeff]*Prod[sys[[2]],Evaluate[b]]];
						rules[[i]] = ReleaseHold[With[{x ={-Evaluate[coeff]*Prod[],i,Prod[Evaluate[b]]},y=q},Alternatives[Pattern[Evaluate[coeff],BlankNullSequence[]]*p,p]:> Hold[Sow[x];y]]],
					True,
						p = Prod[Pattern[Evaluate[a],BlankNullSequence[]],sys[[1]],Pattern[Evaluate[b],BlankNullSequence[]]];
						q = Expand[Evaluate[coeff]*Prod[Evaluate[a],sys[[2]],Evaluate[b]]];
						rules[[i]] = ReleaseHold[With[{x ={-Evaluate[coeff]*Prod[Evaluate[a]],i,Prod[Evaluate[b]]},y=q},Alternatives[Pattern[Evaluate[coeff],BlankNullSequence[]]*p,p]:> Hold[Sow[x];y]]];		
				];
				i = 1,
				i++;
			];
		];
	];
	{DeleteCases[G,Null],Map[ExpandLeft[CollectLeft[#]]&,DeleteCases[cofactors,Null]]}
]


interreduceHead[ideal_,OptionsPattern[{Rules->None,OneSided->False}]]:=
Module[{G,rules,i,s,gi,cofactors,r,lt,sys,a,b,coeff,p,q,newPart,lc,tail,oneSided},
	oneSided = ToLowerCase[OptionValue[OneSided]];
	(*set everything up*)
	G = If[FreeQ[ideal,Prod],
		ToProd/@ideal,
		ideal
	];
	lc = MakeMonic[G];
	cofactors = MapIndexed[{{1/#1*Prod[],#2[[1]],Prod[]}}&,lc];
	
	sys = CreateRedSys[ideal];
	rules = If[OptionValue[Rules] === None,
		Which[oneSided === "left",
				ExtractLeftRules[sys],
			oneSided === "right",
				ExtractRightRules[sys],
			True,
				ExtractRules[sys]
		],
		OptionValue[Rules]
	];
	i = 1; s = Length[G];
	
	(*actual interreduction*)
	While[i <= s,
		If[G[[i]]===Null,i++;Continue[]];
		gi = G[[i]];
		lt = LeadingTermIntern[gi];
		tail = Remainder[gi,lt];
		r = Reap[lt = ToProd[lt[[2]]]/.Drop[rules,{i}]];
		gi = lt + tail;
		If[gi===0,
			rules[[i]] = Null->Null;
			G[[i]] = Null;
			cofactors[[i]] = Null;
			i++,
			If[gi =!= G[[i]],
				lt = LeadingTermIntern[gi];
				(* if rules was used for tracing, adapt cofactors *)
				If[Length[r[[2]]] > 0,
					r = r[[2,1]];
					newPart = Flatten[Table[Map[{Prod[r[[i,1]],#[[1]]],#[[2]],Prod[#[[3]],r[[i,3]]]}&,cofactors[[r[[i,2]]]]],{i,Length[r]}],1];
					cofactors[[i]] = Join[cofactors[[i]],newPart];
				];
				(* make poly monic if it is not already *)
				If[lt[[1]] =!= 1, 
					gi = Expand[1/lt[[1]]*gi]; 
					cofactors[[i]] = (ReplacePart[#,1 -> 1/lt[[1]]*#[[1]]]&/@ cofactors[[i]]);
				];
				G[[i]]=gi;
				sys = CreateRedSys[gi];
				Which[oneSided === "left",
						p = Prod[Pattern[Evaluate[a],BlankNullSequence[]],sys[[1]]];
						q = Expand[Evaluate[coeff]*Prod[Evaluate[a],sys[[2]]]];
						rules[[i]]= ReleaseHold[With[{x ={-Evaluate[coeff]*Prod[Evaluate[a]],i,Prod[]},y=q},Alternatives[Pattern[Evaluate[coeff],BlankNullSequence[]]*p,p]:> Hold[Sow[x];y]]],
					oneSided === "right",
						p = Prod[sys[[1]],Pattern[Evaluate[b],BlankNullSequence[]]];
						q = Expand[Evaluate[coeff]*Prod[sys[[2]],Evaluate[b]]];
						rules[[i]] = ReleaseHold[With[{x ={-Evaluate[coeff]*Prod[],i,Prod[Evaluate[b]]},y=q},Alternatives[Pattern[Evaluate[coeff],BlankNullSequence[]]*p,p]:> Hold[Sow[x];y]]],
					True,
						p = Prod[Pattern[Evaluate[a],BlankNullSequence[]],sys[[1]],Pattern[Evaluate[b],BlankNullSequence[]]];
						q = Expand[Evaluate[coeff]*Prod[Evaluate[a],sys[[2]],Evaluate[b]]];
						rules[[i]] = ReleaseHold[With[{x ={-Evaluate[coeff]*Prod[Evaluate[a]],i,Prod[Evaluate[b]]},y=q},Alternatives[Pattern[Evaluate[coeff],BlankNullSequence[]]*p,p]:> Hold[Sow[x];y]]];		
				];
				i = 1,
				i++;
			];
		];
	];
	{DeleteCases[G,Null],Map[ExpandLeft[CollectLeft[#]]&,DeleteCases[cofactors,Null]]}
]


(* ::Subsection::Closed:: *)
(*Homogeneous & Monomial part*)


SetAttributes[Hom,HoldFirst];

Hom[cofactors_,ideal_,maxiter:_?IntegerQ:10,A:_?MatrixQ:{},OptionsPattern[{Criterion->True,Ignore->0,MaxDeg->Infinity,Parallel->True,Sorted->True}]]:=
Module[{vars,deg,m,n,I,T,TInv,t,tinv,oldOrdering,commutators,J,G,pos,rules,p},
	
	I = ToProd/@ideal;
	vars = WordOrder;
	
	(* If deg \[Equal] {}, then set deg = I_n *)
	If[A === {},
		deg = IdentityMatrix[Length[vars]],
		deg = A
	];
	{m,n} = Dimensions[deg];

	T = Table[Subscript[t, i],{i,m}];
	TInv = Table[Subscript[tinv, i],{i,m}];
	
	(* set up new ordering *)
	If[SortedQ === DegLex,
		oldOrdering = {DegLex,WordOrder};
		SetUpRing[WordOrder,T,TInv],
		oldOrdering = {MultiLex,varSets};
		SetUpRing[Sequence@@(varSets//Reverse),Join[T,TInv]]
	];
	
	(* set up ideal *)
	J = Join[Outer[#1**#2-#2**#1 &,vars,T],Map[Prod[#]-Prod[#//Reverse]&,Subsets[T,{2}]],Outer[#1**#2-#2**#1 &,vars,TInv],Map[Prod[#]-Prod[#//Reverse]&,Subsets[TInv,{2}]],Outer[#1**#2-#2**#1 &,T,TInv],MapThread[1-#1**#2 &,{T,TInv}]]//Flatten;
	(* multiply variables by t's according to deg *)
	rules = {};
	Do[
		p = Table[If[deg[[i,j]] >= 0,
					ConstantArray[T[[j]],deg[[i,j]]],
					ConstantArray[TInv[[j]],-deg[[i,j]]]],
			{j,n}];
		AppendTo[rules,vars[[i]] -> Prod[p,vars[[i]]]];
	,{i,m}];
	I = (I/.rules)//.ExtractRules[CreateRedSys[J]];
	
	(* compute Gr\[ODoubleDot]bner basis and intersect with K<X> *)
	G = Groebner[cofactors,ToNonCommutativeMultiply/@Join[J,I],maxiter,Criterion->OptionValue[Criterion],Ignore->OptionValue[Ignore],MaxDeg->OptionValue[MaxDeg],Parallel->OptionValue[Parallel],Sorted->OptionValue[Sorted]];
	pos = Position[G,g_/;FreeQ[g,Alternatives@@Join[T,TInv]],{1},Heads->False];
	G = Extract[G,pos];
	
	(* process cofactors *)
	cofactors = Map[ToProd,Extract[cofactors,pos],{3}];
	cofactors = DeleteCases[cofactors/.Map[#->1&,Join[T,TInv]],{__,0,__},{2}];
	cofactors = Map[ToNonCommutativeMultiply,cofactors,{3}];
	
	(* Change ordering back *) 
	If[oldOrdering[[1]] === DegLex,
		SetUpRing[oldOrdering[[2]]],
		SetUpRing[Sequence@@(oldOrdering[[2]]//Reverse)]
	];
	
	G
]


SetAttributes[Mon,HoldFirst];

Mon[cofactors_,ideal_,maxiter:_?IntegerQ:10, OptionsPattern[{OneSided->"right",Criterion->True,Ignore->0,MaxDeg->Infinity,Parallel->True,Sorted->True}]]:=
Module[{oneSided,oldOrdering,I,vars,T,TInv,t,tinv,tag,J,G,pos},
	oneSided = ToLowerCase[OptionValue[OneSided]];

	I = ToProd/@ideal;
	vars =  DeleteDuplicates[I//.{Plus -> List, Times[_,v_] :> v, Prod -> List}//Flatten];
	T = Table[Subscript[t, i],{i,Length[vars]}];
	TInv = Table[Subscript[tinv, i],{i,Length[vars]}];
	
	(* set up new ordering *)
	If[SortedQ === DegLex,
		oldOrdering = {DegLex,WordOrder};
		SetUpRing[WordOrder,{tag},Join[T,TInv]],
		oldOrdering = {MultiLex,varSets};
		SetUpRing[Sequence@@(varSets//Reverse),{tag},Join[T,TInv]]
	];
	
	(* set up ideal to factor out and add t variables and tag variable *)
	Switch[oneSided,
			"right",
				J = Join[Outer[#1**#2-#2**#1&,vars,T],MapThread[1-#1**#2&,{T,TInv}]]//Flatten;
				I = I/.Prod[m__]:>Prod[m,Extract[T,Position[vars,#,{1}]&/@{m}]];
				I = Map[Prod[tag,#]&,I],
			"left",
				J = Join[Outer[#1**#2-#2**#1&,vars,TInv],MapThread[1-#1**#2&,{TInv,T}]]//Flatten;
				I = I/.Prod[m__]:>Prod[Extract[T,Position[vars,#,{1}]&/@{m}],m];
				I = Map[Prod[#,tag]&,I],
			_,
				Print["Optional argument OneSided has invalid value: ", oneSided];
				Return[]
	];
	
	(* compute Gr\[ODoubleDot]bner basis and intersect with K<X> *)
	G = Groebner[cofactors,ToNonCommutativeMultiply/@Join[J,I],maxiter,Criterion->OptionValue[Criterion],Ignore->Length[J],MaxDeg->OptionValue[MaxDeg],Parallel->OptionValue[Parallel],Sorted->OptionValue[Sorted]];
	pos = Position[G,g_/;FreeQ[g,Alternatives@@Join[T,TInv]],{1},Heads->False];
	G = Extract[G,pos];
	
	(* process cofactors *)
	cofactors = Map[ToProd,Extract[cofactors,pos],{3}];
	cofactors = DeleteCases[cofactors/.Map[#->1&,Join[T,TInv]],{__,0,__},{2}];
	
	(* Remove the tag variable *)
	G = G/.{NonCommutativeMultiply[l___,tag,r___]:>NonCommutativeMultiply[l,r],tag->1};
	cofactors = cofactors/.{tag->1};
	cofactors = Map[ToNonCommutativeMultiply,cofactors,{3}];
	
	(* Change ordering back *) 
	If[oldOrdering[[1]] === DegLex,
		SetUpRing[oldOrdering[[2]]],
		SetUpRing[Sequence@@(oldOrdering[[2]]//Reverse)]
	];
	
	G
]


(* ::Section::Closed:: *)
(*Module GB computations*)


(* ::Subsection::Closed:: *)
(*Module basis*)


SetModuleBasis[B_List] := (ModuleBasis := B);


(* ::Subsection::Closed:: *)
(*Module Additional stuff*)


SetAttributes[ModuleMakeMonic,HoldFirst];

ModuleMakeMonic[ideal_]:= Module[{lc},
	lc = (ModuleLeadingTermIntern/@ideal)[[All,1]];
	ideal = Expand[ideal/lc];
	lc
]


ModuleCreateRedSys[polies_List]:=
	DeleteCases[Map[ModuleCreateRedSys,polies],{}]


ModuleCreateRedSys[p_]:= Module[{lt,poly},
	poly = ToProd[p];
	If[poly===0,
		Return[{}]
	];
	lt = ModuleLeadingTermIntern[poly];
	{lt[[2]], Expand[-1/lt[[1]]*Remainder[poly,lt]]}
]


(* ::Subsection::Closed:: *)
(*Module Ordering*)


ModuleSortedQ[a_List,b_List]:=Module[{ei,ej,u1,u2,v1,v2,i,j},
	{u1,ei,u2}=a;
	{v1,ej,v2}=b;
	i = Position[ModuleBasis,ei,{1}][[1,1]];
	j = Position[ModuleBasis,ej,{1}][[1,1]];
	If[i < j ,
		True,
		If[i > j,
			False,
			(* alar < blbr iff !(blbr \[LessEqual] alar)*)
			If[!SortedQ[Join[v1,v2],Join[u1,u2]],
				True,
				If[Join[v1,v2]===Join[u1,u2],
					SortedQ[u1,v1],
					False
				]
			]
		]
	]
]


ModuleLeadingTerm[poly_] := ModuleLeadingTermIntern[ToProd[poly]]


ModuleLeadingTermIntern[poly_Prod]:=Module[{ei,lm,m,pos},
	m = List@@poly;
	ei = Intersection[ModuleBasis,m][[1]];
	pos = Position[m,ei,{1}][[1,1]];
	{1,{m[[;;pos-1]],ei,m[[pos+1;;]]}}
]


ModuleLeadingTermIntern[poly_]:=
	If[Head[Expand[poly]] === Times,
		{First[poly],ModuleLeadingTermIntern[Last[poly]][[2]]},
		ModuleFindMaxTerm[Expand[poly]]
	]


ModuleFindMaxTerm[poly_]:= Module[{terms},
	terms = MonomialList[poly];
	Last[Sort[ModuleLeadingTermIntern/@terms,ModuleSortedQ[#1[[2]],#2[[2]]]&]]
]


(* ::Subsection::Closed:: *)
(*Module Ambiguities*)


ModuleAmbiguity[{{v1_,ei_,v2_},i_Integer},{{w1_,ej_,w2_},j_Integer}]:=
Module[{min1,min2,max1,max2,amb},
	amb = {};
	If[ei === ej,
		min1 = Min[Length[v1],Length[w1]];
		min2 = Min[Length[v2],Length[w2]];
		If[v1[[-min1;;]] === w1[[-min1;;]] && v2[[;;min2]] === w2[[;;min2]],
			max1 = MaximalBy[{v1,w1},Length,1][[1]];
			max2 = MaximalBy[{v2,w2},Length,1][[1]];
			Which[
				{max1,max2}==={v1,v2},
					amb = ModuleAmbiguity[{max1,ei,max2},{},{},v1[[;;-min1-1]],v2[[min2+1;;]],{i,j}],
				{max1,max2}==={v1,w2},
					amb = ModuleAmbiguity[{max1,ei,max2},{},w2[[min2+1;;]],v1[[;;-min1-1]],{},{i,j}],
				{max1,max2}==={w1,v2},
					amb = ModuleAmbiguity[{max1,ei,max2},w1[[;;-min1-1]],{},{},v2[[min2+1;;]],{i,j}],
				{max1,max2}==={w1,w2},
					amb = ModuleAmbiguity[{max1,ei,max2},w1[[;;-min1-1]],w2[[min2+1;;]],{},{},{i,j}]
			];
		];
	];
	amb
]


ModuleGenerateAmbiguities[l_List,newpart_List,maxdeg_:Infinity,OptionsPattern[{Parallel->True}]]:= 
	Select[
		If[OptionValue[Parallel],
			Join[
				Parallelize[Outer[If[#1[[2]] < #2[[2]],ModuleAmbiguity[#1,#2],{}]&,l,newpart,1],DistributedContexts->Automatic],
				Parallelize[Outer[If[#1[[2]] < #2[[2]],ModuleAmbiguity[#1,#2],{}]&,newpart,newpart,1],DistributedContexts->Automatic]
			],
			Print["here"];
			Join[
				Outer[If[#1[[2]] < #2[[2]],ModuleAmbiguity[#1,#2],{}]&,l,newpart,1],
				Outer[If[#1[[2]] < #2[[2]],ModuleAmbiguity[#1,#2],{}]&,newpart,newpart,1]
			]
		]//Flatten, Length[#[[1]]] <= maxdeg &]


ModuleGenerateAmbiguities[l_List,maxdeg_:Infinity,OptionsPattern[{Parallel->True}]] := 
	Select[
		If[OptionValue[Parallel],
			Parallelize[Outer[If[#1[[2]] < #2[[2]],ModuleAmbiguity[#1,#2],{}]&,l,l,1],DistributedContexts->Automatic],
			Outer[If[#1[[2]] < #2[[2]],ModuleAmbiguity[#1,#2],{}]&,l,l,1]
		]//Flatten, Length[#[[1]]] <= maxdeg &]


(* ::Subsection::Closed:: *)
(*Module S-Polynomials*)


ModuleSPoly[amb:_ModuleAmbiguity,fi_,fj_]:=
Module[{v1,v2,w1,w2},
	{v1,v2,w1,w2} = Map[Prod@@#&,List@@amb[[2;;-2]]];
	{Prod[v1,fi[[2]],v2] - Prod[w1,fj[[2]],w2],
	{{-v1,amb[[6,1]],v2},{w1,amb[[6,2]],w2}}}
]


(* ::Subsection::Closed:: *)
(*Module Groebner basis*)


(* ::Text:: *)
(*Implementation of the Buchberger algorithm to compute a (partial) Groebner basis of the module 'M' with at most 'maxiter' iterations being executed (default: 10). The return value is a  set of module elements G = {f1,...,fn,g1,...gm} consisting of the elements f1,....,fn from 'M' and new elements g1,...,gm. For every element g\in G a list forming a linear combination of g consisting of elements from 'M' and certain cofactors is saved in the list cofactors.*)
(**)
(*computation.*)
(*	- Ignore (default: 0): A non-negative integer that determines how many elements of the input will be ignored during the first computation of the ambiguities. *)
(*	- MaxDeg (default: Infinity): Only ambiguities with degree smaller than or equal to MaxDeg will be considered during the Groebner basis computation (larger ambiguities are simply ignored). *)
(*	- Parallel (default: True): Determines whether the computations for which it is possible, are executed in parallel (which speeds up the computation) or in series.*)
(*	- Sorted (default: True):  Sorts the ambiguities before processing in ascending order. This speeds up the computation but results in a different (partial) Groebner basis.*)
(*	- IterCount (default: 0): defines from which number the iterations are counted (only relevant for the printed information)*)
(**)


SetAttributes[ModuleGroebner,HoldFirst]

ModuleGroebner[cofactors_,M_, maxiter:_?IntegerQ:10, OptionsPattern[{Ignore->0,MaxDeg->Infinity,Parallel->True,Sorted->True,IterCount->0}]]:=
Module[{lc,count,spol,lt,p,h,G,r,t1,t2,rules,sorted,oldlength,parallel,hrule,maxdeg,intern,i,itercount},
	sorted = OptionValue[Sorted];
	parallel = OptionValue[Parallel];
	maxdeg = OptionValue[MaxDeg];
	itercount = OptionValue[IterCount];
	intern = !FreeQ[ideal,Prod]; 
	
	If[Head[cofactors] =!= List,
		cofactors = {};
	];

	If[intern,
		G = M,
		G = ToProd/@M;
		lc = ModuleMakeMonic[G];
		cofactors = MapIndexed[{{1/#1*Prod[],#2[[1]],Prod[]}}&,lc];
	];
	
	G = ModuleCreateRedSys[G];
	oldlength = OptionValue[Ignore];
	t1 = 0; t2 = 0; count = 0;
	If[VerboseOperatorGB > 0,
		Print["G has ", Length[G]," elements in the beginning.\n"]
	];

	rules = ExtractRules[G];
	spol = {0};

	While[Length[spol] > 0 && count < maxiter,
		count++;
		If[VerboseOperatorGB > 0,
			Print["Starting iteration ", count + itercount ,"..."]
		];
		spol = ModuleCheckResolvability[G,oldlength,MaxDeg->maxdeg,Sorted->sorted,Parallel->parallel];
		oldlength = Length[G];
		
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
				lt = ModuleLeadingTermIntern[h];
				If[lt[[1]] =!= 1, 
					h = Expand[1/lt[[1]]*h]; p[[2]] = (ReplacePart[#,1 -> 1/lt[[1]]*#[[1]]]&/@ p[[2]])
				];
				hrule = ModuleCreateRedSys[h];
				AppendTo[G,hrule];
				AppendTo[cofactors,p[[2]]]; 
				AppendTo[rules,ExtractRule[hrule,Length[G]]];
			];
			i--;
		,{p,spol}];,i];][[1]];
		
		If[VerboseOperatorGB > 1, 
			Print["The second reduction took ", t1]
		];
		If[VerboseOperatorGB > 0,
			Print["Iteration ",count + itercount, " finished. G has now ", Length[G]," elements\n"]
		];
	];
	
	G = ToPoly[G];
	
	If[intern,
		G,
		RewriteGroebner[cofactors,M];
		ToNonCommutativeMultiply[G]
	]
]


(* ::Text:: *)
(*ModuleCheckResolvability returns all S-polynomials from the reduction system sys which can not be reduced to zero. Additionally, for each S-polynomial a list containing the linear combination how the S-polynomial was generated from the elements of sys is returned.*)
(*For a description of the OptionPatterns see the documentation of the ModuleGroebner method.*)


ModuleCheckResolvability[sys_,oldlength:_?IntegerQ:0,OptionsPattern[{MaxDeg->Infinity,Parallel->True,Sorted->True}]]:=
Module[{amb,spol,rules,parallel,words,r,t1,t2,t3,str},
	parallel = OptionValue[Parallel];

	(*generate ambiguities*)
	words = ExtractReducibleWords[sys];
	t1 = AbsoluteTiming[
		amb = ModuleGenerateAmbiguities[words[[;;oldlength]],words[[oldlength+1;;]],OptionValue[MaxDeg],Parallel->parallel];
	][[1]];
	If[VerboseOperatorGB > 0,
		str = ToString[Length[amb]] <> " ambiguities in total" <> If[VerboseOperatorGB > 1, " (computation took " <> ToString[t1] <> ")",""];
		Print[str];
	];
	
	(*process ambiguities*)
	If[OptionValue[Sorted],
		amb = SortBy[amb,Length[#[[1]]]&];
	];
	
	(*generate S-polynomials*)
	t2 = AbsoluteTiming[
	spol = DeleteCases[
		If[parallel,
			ParallelMap[ModuleSPoly[#,sys[[#[[6,1]]]],sys[[#[[6,2]]]]]&,amb,DistributedContexts->Automatic,Method->"CoarsestGrained"],
			Map[ModuleSPoly[#,sys[[#[[6,1]]]],sys[[#[[6,2]]]]]&,amb]
		]
	,{0,_}];
	][[1]];
	
	If[VerboseOperatorGB > 1,
		Print["Generating S-polys: ",t2 ," (",Length[spol], " in total)"]
	];
	
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
	If[VerboseOperatorGB > 1,
		Print["Reducing S-polys: ",t3, " (",Length[spol], " remaining)"]
	];
	If[OptionValue[Sorted],
		If[$VersionNumber < 12,
			Sort[spol,ModuleSortedQ[ModuleLeadingTermIntern[#1[[1]]][[2]],ModuleLeadingTermIntern[#2[[1]]][[2]]]&],
			SortBy[spol,ModuleLeadingTermIntern[#[[1]]][[2]]&,ModuleSortedQ]
		],
		spol
	]
]


(* ::Section::Closed:: *)
(*Quiver*)


(* ::Subsection::Closed:: *)
(*Quiver data structure & basic functionality *)


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


(* ::Subsection::Closed:: *)
(*Constructing quivers from polynomials*)


TrivialQuiver[F_]:= Module[{vars,v},
	vars =  DeleteDuplicates[(ToProd/@F)//.{Plus -> List, Times[_,v_] :> v, Prod -> List}//Flatten];
	Table[{v,1,1},{v,vars}]
]


QuiverFromPolynomials[FF_]:=
Module[{F,vars,monomials,m,f,i,rowIdx,sourceIdx,sourceIdx2 ,targetIdx,targetIdx2,M},
	F = ToProd/@FF;
	vars = DeleteDuplicates[F//.{Plus->List,Times[_,m_]:>m,Prod->List}//Flatten];
	vars = Thread[Rule[vars,Range[Length[vars]]]];
	monomials = Monomials/@F;

	(* special case: a polynomial has a constant term *)
	If[MemberQ[monomials,{},Infinity],
		Return[Map[{#,1,1}&,vars[[All,1]]]]
	];

	(* usual case: no constant term *)
	(* make sure each individual monomial is compatible *)
	rowIdx = 1;
	M = Reap[
	Do[
		Do[
			sourceIdx = 2*( m[[i]]/.vars)-1;
			targetIdx =2* ( m[[i+1]]/.vars);
			Sow[{{rowIdx,sourceIdx}->1,{rowIdx,targetIdx}->-1}];
			rowIdx++;
		,{i,Length[m]-1}]
	,{m,Flatten[monomials,1]}];

	(* make sure the monomials of each polynomial fit together *)
	Do[
		Do[
			sourceIdx = 2*( f[[i,-1]]/.vars)-1;
			targetIdx = 2* ( f[[i,1]]/.vars);
			sourceIdx2 = 2*( f[[i+1,-1]]/.vars)-1;
			targetIdx2 = 2* ( f[[i+1,1]]/.vars);
			If[sourceIdx =!= sourceIdx2,
				Sow[{{rowIdx,sourceIdx}->1,{rowIdx,sourceIdx2}->-1}];
				rowIdx++;
			];
			If[targetIdx =!= targetIdx2,
				Sow[{{rowIdx,targetIdx}->1,{rowIdx,targetIdx2}->-1}];
				rowIdx++;
			];
		,{i,Length[f]-1}]
	,{f,monomials}];
	][[2]];
	If[Length[M]===0,
		Print["special case - not implemented yet"];
		Return[]
	];
	M = NullSpace[SparseArray[Flatten[M,2],{rowIdx-1,2*Length[vars]}]]//Abs;
	Table[{vars[[i,1]],Position[M[[All,2*i-1]],1][[1,1]],Position[M[[All,2*i]],1 ][[1,1]]},{i,Length[vars]}]
];


(* ::Subsection::Closed:: *)
(*Compatibility checks*)


(* ::Text:: *)
(*Returns the signature of a polynomial w.r.t. to a quiver Q (not necessarily with unique labels)*)


QSignature[l_List,Q:Quiver]:= Map[QSignature[#,Q]&,l]


QSignature[p_,Q:Quiver]:=
Module[{monomials,m,results,i,begin,end,comb,vertices,x},
	monomials = Monomials[ToProd[p]];
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


CompatibleQ[p_,Q:Quiver] := QSignature[p,Q]=!={}


UniformlyCompatibleQ[p_,Q:Quiver]:= 
Module[{monomials,signatures},
	If[!CompatibleQ[p,Q],
		Return[False]
	];
	monomials = MonomialList[p];
	signatures = Map[QSignature[#,Q]&,monomials];
	AllTrue[signatures,#===signatures[[1]]&]
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

QCompletion[cofactors_,ideal_,Q:Quiver,maxiter:_?IntegerQ:10, OptionsPattern[{Criterion->True,Ignore->0,MaxDeg->Infinity,Parallel->True,Sorted->True,IterCount->0}]]:=
Module[{lc,count,spol,lt,p,h,G,r,t1,t2,rules,sorted,oldlength,parallel,hrule,maxdeg,intern,criterion,i,itercount},
	sorted = OptionValue[Sorted];
	parallel = OptionValue[Parallel];
	maxdeg = OptionValue[MaxDeg];
	criterion = OptionValue[Criterion];
	itercount = OptionValue[IterCount];
	intern = !FreeQ[ideal,Prod]; 

	If[intern,
		G = ideal,
		G = ToProd/@ideal;
		lc = MakeMonic[G];
		cofactors = MapIndexed[{{1/#1*Prod[],#2[[1]],Prod[]}}&,lc];
	];
	
	G = CreateRedSys[G];
	oldlength = OptionValue[Ignore];
	t1 = 0; t2 = 0; count = 0;
	If[VerboseOperatorGB > 0,
		Print["G has ", Length[G]," elements in the beginning.\n"]
	];

	rules = ExtractRules[G];
	spol = {0};

	While[Length[spol] > 0 && count < maxiter,
		count++;
		If[VerboseOperatorGB > 0,
			Print["Starting iteration ", count + itercount ,"..."]
		];
		spol = CheckQResolvability[G,Q,oldlength,Criterion->criterion,MaxDeg->maxdeg,Sorted->sorted,Parallel->parallel];
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
		
		If[VerboseOperatorGB > 1, 
			Print["The second reduction took ", t1];
		];
		If[VerboseOperatorGB > 0,
			Print["Iteration ",count + itercount, " finished. G has now ", Length[G]," elements\n"]
		];
	];
	
	G = ToPoly[G];
	
	If[intern,
		G,
		RewriteGroebner[cofactors,ideal];
		ToNonCommutativeMultiply[G]
	]
]


CheckQResolvability[sys_, Q:Quiver, oldlength:_?IntegerQ:0,OptionsPattern[{Criterion->True,MaxDeg->Infinity,Parallel->True,Sorted->True}]]:=
Module[{amb,spol,rules,parallel,words,r,t1,t2,t3,str},
	parallel = OptionValue[Parallel];

	(*generate ambiguities*)
	words = ExtractReducibleWords[sys];
	t1 = AbsoluteTiming[
		amb = GenerateAmbiguities[words[[;;oldlength]],words[[oldlength+1;;]],OptionValue[MaxDeg],Parallel->parallel];
	][[1]];
	
	(* remove ambiguities whose source is not compatible *)
	amb = Select[amb,CompatibleQ[Prod@@#[[1]],Q]&];
	
	If[VerboseOperatorGB > 0,
		str = ToString[Length[amb]] <> " ambiguities in total" <> If[VerboseOperatorGB > 1, " (computation took " <> ToString[t1] <> ")", ""];
		Print[str];
	];
	
	(*process ambiguities*)
	If[OptionValue[Criterion],
		amb = DeleteRedundant[amb,words];
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
	If[VerboseOperatorGB > 1,
		Print["Generating S-polys: ",t2 ," (",Length[spol], " in total)"]
	];
	
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
	If[VerboseOperatorGB > 1,
		Print["Reducing S-polys: ",t3, " (",Length[spol], " remaining)"]
	];
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


(* ::Section::Closed:: *)
(*Certifying operator identities*)


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
(*Cancellability*)


ApplyLeftCancellability[II_,aabb_,aa_,OptionsPattern[{MaxIter->1, Algorithm->"subalgebra", Vars->{}}]]:=
Module[{G,f,I,ab,cancel,rule,a,b,c,vars,subalg,Q,abx,pos},
	If[FreeQ[II,Prod] ||\[NonBreakingSpace]FreeQ[aabb,Prod] || FreeQ[aa,Prod],
		I = ToProd/@II; ab = ToProd[aabb]; a = ToProd[aa],
		I = II; ab = aabb; a = aa;
	];
	(* all variables that should be considered during the computation *)
	If[Length[OptionValue[Vars]] > 0,
		vars = OptionValue[Vars],
		vars = WordOrder
	];
	(* compute Groebner basis of intersection either of ideals or of ideal with subalgebra *)
	Switch[OptionValue[Algorithm],
		"two-sided",
		G = Intersect[I,{ab},MaxIter->OptionValue[MaxIter],Vars->vars],
		"one-sided",
		Q = QuiverFromPolynomials[Join[I,{ab}]];
		G = IntersectRightIdeal[I,{ab},Q,Length->OptionValue[MaxIter],Vars->vars],
		_,
		subalg = Join[{ab},ToProd/@vars];
		G = IntersectSubalgebra[I,subalg,MaxIter->OptionValue[MaxIter]];
	];
	G = ToProd/@G;
	(* only pick those elements which are in the intersection when (ab) is considered as a right ideal *)
	rule = ExtractRightRules[CreateRedSys[{ab}]];
	pos = Position[G//.rule,0];
	abx = Extract[G,pos];
	
	(* transform abx into bx *)
	rule = ExtractRightRules[CreateRedSys[{a}]];
	cancel = Cases[Map[Reap[#//.rule]&,abx],{0,{c__}}:>c];
	cancel = cancel /.{{Prod[],_,c_} :> -c, {-Prod[],_,c_} :> c};
	ToNonCommutativeMultiply/@(Plus@@#&/@cancel)
]


ApplyRightCancellability[II_,aabb_,bb_,OptionsPattern[{MaxIter->1, Algorithm->"subalgebra", Vars->{}}]]:=
Module[{G,f,I,ab,cancel,rule,pos,a,b,c,xab,vars,subalg,Q},
	If[FreeQ[II,Prod] ||\[NonBreakingSpace]FreeQ[aabb,Prod] || FreeQ[aa,Prod],
		I = ToProd/@II; ab = ToProd[aabb]; b = ToProd[bb],
		I = II; ab = aabb; b = bb;
	];
	(* all variables that should be considered during the computation *)
	If[Length[OptionValue[Vars]] > 0,
		vars = OptionValue[Vars],
		vars = WordOrder
	];
	(* compute Groebner basis of intersection either of ideals or of ideal with subalgebra *)
	Switch[OptionValue[Algorithm],
		"two-sided",
		G = Intersect[I,{ab},MaxIter->OptionValue[MaxIter],Vars->vars],
		"one-sided",
		Q = QuiverFromPolynomials[Join[I,{ab}]];
		G = IntersectLeftIdeal[I,{ab},Q,Length->OptionValue[MaxIter],Vars->vars],
		_,
		subalg = Join[{ab},ToProd/@vars];
		G = IntersectSubalgebra[I,subalg,MaxIter->OptionValue[MaxIter]];
	];
	G = ToProd/@G;
	(* only pick those elements which are in the intersection when (ab) is considered as a left ideal *)
	rule = ExtractLeftRules[CreateRedSys[{ab}]];
	pos = Position[G//.rule,0];
	xab = Extract[G,pos];
	(* transform xab into xa *)
	rule = ExtractLeftRules[CreateRedSys[{b}]];
	cancel = Cases[Map[Reap[#//.rule]&,xab],{0,{c__}}:>c];
	cancel = cancel /.{c_,_,Prod[]} :> -c;
	ToNonCommutativeMultiply/@(Plus@@#&/@cancel)
]


ExtractRightRules[sys_]:=
Module[{b,coeff,i,p,q,f},
	b=Unique[];coeff=Unique[];
	MapIndexed[(i = #2[[1]]; f = #1;
		p = Prod[f[[1]],Pattern[Evaluate[b],BlankNullSequence[]]];
		q = Expand[Evaluate[coeff]*Prod[f[[2]],Evaluate[b]]];
		With[{x ={-Evaluate[coeff]*Prod[],i,Prod[Evaluate[b]]},y=q},Alternatives[Pattern[Evaluate[coeff],BlankNullSequence[]]*p,p]:> Hold[Sow[x];y]])&
	,sys]//ReleaseHold
]


ExtractLeftRules[sys_]:=
Module[{a,coeff,i,p,q,f},
	a=Unique[];coeff=Unique[];
	MapIndexed[(i = #2[[1]]; f = #1;
		p = Prod[Pattern[Evaluate[a],BlankNullSequence[]],f[[1]]];
		q = Expand[Evaluate[coeff]*Prod[Evaluate[a],f[[2]]]];
		With[{x ={-Evaluate[coeff]*Prod[Evaluate[a]],i,Prod[]},y=q},Alternatives[Pattern[Evaluate[coeff],BlankNullSequence[]]*p,p]:> Hold[Sow[x];y]])&
	,sys]//ReleaseHold
]


(* ::Subsection::Closed:: *)
(*Certify*)


Certify[assumptionsInput_List,claimsInput_,Q:Quiver,OptionsPattern[{MaxIter->10,MaxDeg->Infinity,MultiLex->False,Parallel->True,Sorted->True,Criterion->True}]]:=
 Module[{maxiter,N,shuffle,normalForm,cofactors,cofactorsReduction,G,sigAssump,sigClaim,certificate,toIgnore,toIgnoreOld,alreadyReduced,i,j,knowns,unknowns,t,assumptions,claims,count},
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
	If[VerboseOperatorGB > 0,
		Print["Using the following monomial ordering:"]
		VerboseOperatorGB++;
	];
	If[OptionValue[MultiLex],
		knowns = DeleteDuplicates[Cases[Q[[All,1]],Alternatives@@Flatten[Map[#/.Alternatives[Plus,Times,NonCommutativeMultiply]->List&,claims]]]];
		unknowns = DeleteDuplicates[Cases[Q[[All,1]],var_/;!MemberQ[knowns,var]]];
		SetUpRing[knowns,unknowns],
		SetUpRing[DeleteDuplicates[Q[[All,1]]]]
	];
	If[VerboseOperatorGB > 1,
		VerboseOperatorGB--;
	];

	(*interreduce the generators*)
	{G,cofactors} = Interreduce[assumptions];
	N = Length[G];
	If[VerboseOperatorGB > 0,
		Print["\nInterreduced the input from ", Length[assumptions], " polynomials to ", Length[G], ".\n"]
	];
	
	(*compute the Groebner basis and reduce the claims*)
	If[VerboseOperatorGB > 0,
		Print["Computing a (partial) Groebner basis and reducing the claim...\n"]
	];
	(*do computation iteratively*)
	toIgnore = 0;
	i = 1;
	alreadyReduced = Table[False,{j,Length[claims]}];
	certificate = Table[{},{j,Length[claims]}];
	While[MemberQ[alreadyReduced,False] && i <= maxiter,
		toIgnoreOld = Length[G];
		G = Groebner[cofactors,G,1,Ignore->toIgnore,MaxDeg->OptionValue[MaxDeg],Parallel->OptionValue[Parallel],Sorted->OptionValue[Sorted],Criterion->OptionValue[Criterion],IterCount->i-1];
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
	normalForm = ToNonCommutativeMultiply/@ToProd/@normalForm;
	If[Head[claimsInput] =!= List, 
		sigClaim = sigClaim[[1]]; certificate = certificate[[1]]; normalForm = normalForm[[1]];
	];
	
	(* print info *)
	If[VerboseOperatorGB > -1,
		If[MemberQ[alreadyReduced,False],
			Print["Failed! Not all claims could be reduced to 0."],
			Print["Done! All claims were successfully reduced to 0."];
		];
	];
	
	If[VerboseOperatorGB > 0,
		(* full info *)
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


(* ::Section::Closed:: *)
(*Diagram chasing*)


(* Find fx in (G) and return x *)
ApplyMono[G_,Q_,f_,OptionsPattern[{Algorithm->"one-sided",MaxIter->3}]]:=
Module[{mono,pos,x,g},
	mono = ToProd/@ApplyLeftCancellability[G,f,f,MaxIter->OptionValue[MaxIter],Algorithm->OptionValue[Algorithm]];
	mono = Select[mono,CompatibleQ[Prod[f,#],Q]&];
	pos = Position[mono//.ExtractRules[CreateRedSys[G]],x_/;x=!=0,1,Heads->False];
	mono = Extract[mono,pos];
	If[VerboseOperatorGB > 0,
		Do[
			Print["Added ",ToNonCommutativeMultiply[g], " because ",f, " is mono."];
		,{g,mono}];
	];
	mono
]


ApplyEpic[Q_,f_,OptionsPattern[{AlreadyFound->{}}]]:=
Module[{alreadyFound,codomainF,vars,epic,y,e},
	alreadyFound = OptionValue[AlreadyFound];
	
	(* Find variables v  with same codomain as f *)
	codomainF = Cases[Q,{f,_,codom_}:>codom][[1]];
	vars = Cases[Q,{v_,_,codom_}:>v/;codom === codomainF && v =!= f ];
	vars = Complement[vars,alreadyFound];
	alreadyFound = Join[alreadyFound,vars];

	(* add ve - fy *)
	y = Table[Unique["y"],{i,Length[vars]}];
	e = Table[Unique["e"],{i,Length[vars]}];
	epic = MapThread[#1**#2 - f**#3&,{vars,e,y}];
	If[VerboseOperatorGB > 0,
		Map[Print["Added ",#, " because ",f, " is epic."]&,epic];
	];
	SetUpRing[Sequence@@Reverse[varSets],Join[y,e]];
	{epic,alreadyFound}
]


(* Find elements of the form gx in (G) and return xe = fy *)
ApplyExact[G_,Q_,{f_,g_},OptionsPattern[{Algorithm->"one-sided",MaxIter->3,AlreadyFound->{}}]]:=
Module[{alreadyFound,x,y,e,i},
	alreadyFound = OptionValue[AlreadyFound];
	(* find elements of the form gx = 0 *)
	x = DeleteDuplicates[DeleteCases[ApplyLeftCancellability[G,g,g,MaxIter->OptionValue[MaxIter],Algorithm->OptionValue[Algorithm]],f]];
	x = Complement[x,alreadyFound];
	alreadyFound = Join[alreadyFound,x];
	x = Select[x,CompatibleQ[Prod[g,#],Q]&];
	If[Length[x] > 0,
		(* x can be factored as xe = fy *)
		y =Table[Unique["y"],Length[x]];
		e = Table[Unique["e"],Length[x]];
		If[VerboseOperatorGB > 0,
			Do[
				Print["Factored ", x[[i]], " into ",x[[i]]**e[[i]]," = ",f** y[[i]]," because of exactness of (",f,",",g,")"];
			,{i,Length[x]}];
		];
		SetUpRing[Sequence@@Reverse[varSets],Join[y,e]];
		x = ToProd/@MapThread[#1**#2 -f**#3&,{x,e,y}];
	];
	{x,alreadyFound}
]


DiagramChase[diagram_,iter_:3,OptionsPattern[{ExactAt -> {},Mono->{},Epic->{},Algorithm->"one-sided",MaxIter->3}]]:=
Module[{G,Q,exactAt,mono,epic,algorithm,maxIter,cofactors,alreadyFoundEpic,alreadyFoundExact,e,i,newEpic,newExact},
	exactAt = OptionValue[ExactAt];
	mono = OptionValue[Mono];
	epic = OptionValue[Epic];
	algorithm = OptionValue[Algorithm];
	maxIter = OptionValue[MaxIter];
	
	alreadyFoundEpic= {}; alreadyFoundExact= {};

	G = ToProd/@diagram;

	Q = QuiverFromPolynomials[G];

	Do[
		Print["Starting iteration ",i,"..."];

		(* apply mono *)
		G = Join[G,Map[ApplyMono[G,Q,#,Algorithm->algorithm,MaxIter->maxIter]&,mono]//Flatten];

		(* apply epic *)
		Do[
			{newEpic,alreadyFoundEpic} = ApplyEpic[Q,e,AlreadyFound->alreadyFoundEpic];
			G = Join[G,newEpic];
			Q = QuiverFromPolynomials[G];
		,{e,epic}];

		(* apply exactness *)
		Do[
			{newExact,alreadyFoundExact}= ApplyExact[G,Q,e,Algorithm->algorithm,MaxIter->maxIter,AlreadyFound->alreadyFoundExact];
			G = Join[G,newExact];
			Q = QuiverFromPolynomials[G];
		,{e,exactAt}];

		(* Groebner basis computation *)
		G = Groebner[cofactors,G,3];

		(* end with new line *)
		Print[];
	,{i,iter}];

	(* return G *)
	ToNonCommutativeMultiply/@G
]


(* ::Section::Closed:: *)
(*End*)


Copyright[a_String,b___String]:= Print[StringJoin[Prepend[{"\n",#}&/@{b},a]]]


Copyright[
    "Package OperatorGB version 1.3.0",
    "Copyright 2019, Institute for Algebra, JKU",
    "by Clemens Hofstadler, clemens.hofstadler@jku.at"];


End[]


EndPackage[]
