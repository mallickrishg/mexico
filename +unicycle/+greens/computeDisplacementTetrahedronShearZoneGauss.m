function [u1,u2,u3]=computeDisplacementTetrahedronShearZoneGauss( ...
    x1,x2,x3,A,B,C,D,e11,e12,e13,e22,e23,e33,nu,varargin)
% function COMPUTEDISPLACEMENTTETRAHEDRONSHEARZONEGAUSS computes the
% displacement field associated with deforming tetrahedral strain volume
% considering the following geometry using the Gauss-Legendre numerical
% quadrature.
%
%                      / North (x1)
%                     /
%        surface     /
%      -------------+-------------- East (x2)
%                  /|
%                 / |     + A
%                /  |    /  .
%                   |   /     .
%                   |  /        .
%                   | /           .
%                   |/              + B
%                   /            .  |
%                  /|          /    |
%                 / :       .       |
%                /  |    /          |
%               /   : .             |
%              /   /|               |
%             / .   :               |
%            +------|---------------+
%          C        :                 D
%                   |
%                   Down (x3)
%
%
% Input:
% x1, x2, x3         north, east coordinates and depth of the observation point,
% A, B, C, D         north, east, and depth coordinates of the vertices,
% eij                source strain component 11, 12, 13, 22, 23 and 33
%                    in the strain volume,
% nu                 Poisson's ratio in the half space.
%
% Option:
% 'N',integer        the number of Gauss integration points
%
% Output:
% u1                 displacement component in the north direction,
% u2                 displacement component in the east direction,
% u3                 displacement component in the down direction.
%
% Author: Sylvain Barbot (sbarbot@ntu.edu.sg) - Feb 9, 2018, Los Angeles.

assert(min(x3(:))>=0,'depth must be positive.');

% process optional input
p = inputParser;
p.addParameter('n',15,@validateN);
p.parse(varargin{:});
optionStruct = p.Results;

% Lame parameter
lambda=2*nu/(1-2*nu);

% isotropic strain
ekk=e11+e22+e33;

% unit normal vectors
nA=cross(C-B,D-B); % BCD
nB=cross(D-C,A-C); % CDA
nC=cross(A-D,B-D); % DAB
nD=cross(B-A,C-A); % ABC

nA=nA/norm(nA);
nB=nB/norm(nB);
nC=nC/norm(nC);
nD=nD/norm(nD);

% check that unit vectors are pointing outward
if (nA'*(A(:)-(B(:)+C(:)+D(:))/3))>0
    nA=-nA;
end
if (nB'*(B(:)-(C(:)+D(:)+A(:))/3))>0
    nB=-nB;
end
if (nC'*(C(:)-(D(:)+A(:)+B(:))/3))>0
    nC=-nC;
end
if (nD'*(D(:)-(A(:)+B(:)+C(:))/3))>0
    nD=-nD;
end

% area of triangle ABC
ABC=norm(cross(C-A,B-A))/2;
% area of triangle BCD
BCD=norm(cross(D-B,C-B))/2;
% area of triangle CDA
CDA=norm(cross(A-C,D-C))/2;
% area of triangle DAB
DAB=norm(cross(B-D,A-D))/2;

% Radii
r1=@(y1,y2,y3) sqrt((x1-y1).^2+(x2-y2).^2+(x3-y3).^2);
r2=@(y1,y2,y3) sqrt((x1-y1).^2+(x2-y2).^2+(x3+y3).^2);

% Green's functions
G11=@(y1,y2,y3) 1/(16*pi*(1-nu))*( ...
    (3-4*nu)./r1(y1,y2,y3)+1./r2(y1,y2,y3)+(x1-y1).^2./r1(y1,y2,y3).^3 ...
    +(3-4*nu)*(x1-y1).^2./r2(y1,y2,y3).^3+2*x3.*y3.*(r2(y1,y2,y3).^2-3*(x1-y1).^2)./r2(y1,y2,y3).^5 ...
    +4*(1-2*nu)*(1-nu)*(r2(y1,y2,y3).^2-(x1-y1).^2+r2(y1,y2,y3).*(x3+y3))./(r2(y1,y2,y3).*(r2(y1,y2,y3)+x3+y3).^2)...
    );
G12=@(y1,y2,y3) (x1-y1).*(x2-y2)/(16*pi*(1-nu)).*( ...
    1./r1(y1,y2,y3).^3+(3-4*nu)./r2(y1,y2,y3).^3-6*x3.*y3./r2(y1,y2,y3).^5 ...
    -4*(1-2*nu)*(1-nu)./(r2(y1,y2,y3).*(r2(y1,y2,y3)+x3+y3).^2) ...
    );
G13=@(y1,y2,y3) (x1-y1)/(16*pi*(1-nu)).*( ...
    (x3-y3)./r1(y1,y2,y3).^3+(3-4*nu)*(x3-y3)./r2(y1,y2,y3).^3 ...
    -6*x3.*y3.*(x3+y3)./r2(y1,y2,y3).^5+4*(1-2*nu)*(1-nu)./(r2(y1,y2,y3).*(r2(y1,y2,y3)+x3+y3)) ...
    );
G21=@(y1,y2,y3) (x1-y1).*(x2-y2)/(16*pi*(1-nu)).*( ...
    1./r1(y1,y2,y3).^3+(3-4*nu)./r2(y1,y2,y3).^3-6*x3.*y3./r2(y1,y2,y3).^5 ...
    -4*(1-2*nu)*(1-nu)./(r2(y1,y2,y3).*(r2(y1,y2,y3)+x3+y3).^2) ...
    );
G22=@(y1,y2,y3) 1/(16*pi*(1-nu))*( ...
    (3-4*nu)./r1(y1,y2,y3)+1./r2(y1,y2,y3)+(x2-y2).^2./r1(y1,y2,y3).^3 ...
    +(3-4*nu)*(x2-y2).^2./r2(y1,y2,y3).^3+2*x3.*y3.*(r2(y1,y2,y3).^2-3*(x2-y2).^2)./r2(y1,y2,y3).^5 ...
    +4*(1-2*nu)*(1-nu)*(r2(y1,y2,y3).^2-(x2-y2).^2+r2(y1,y2,y3).*(x3+y3))./(r2(y1,y2,y3).*(r2(y1,y2,y3)+x3+y3).^2)...
    );
G23=@(y1,y2,y3) (x2-y2)/(16*pi*(1-nu)).*( ...
    (x3-y3)./r1(y1,y2,y3).^3+(3-4*nu)*(x3-y3)./r2(y1,y2,y3).^3 ...
    -6*x3.*y3.*(x3+y3)./r2(y1,y2,y3).^5+4*(1-2*nu)*(1-nu)./(r2(y1,y2,y3).*(r2(y1,y2,y3)+x3+y3)) ...
    );
G31=@(y1,y2,y3) (x1-y1)/(16*pi*(1-nu)).*( ...
    (x3-y3)./r1(y1,y2,y3).^3+(3-4*nu)*(x3-y3)./r2(y1,y2,y3).^3 ...
    +6*x3.*y3.*(x3+y3)./r2(y1,y2,y3).^5-4*(1-2*nu)*(1-nu)./(r2(y1,y2,y3).*(r2(y1,y2,y3)+x3+y3)) ...
    );
G32=@(y1,y2,y3) (x2-y2)/(16*pi*(1-nu)).*( ...
    (x3-y3)./r1(y1,y2,y3).^3+(3-4*nu)*(x3-y3)./r2(y1,y2,y3).^3 ...
    +6*x3.*y3.*(x3+y3)./r2(y1,y2,y3).^5-4*(1-2*nu)*(1-nu)./(r2(y1,y2,y3).*(r2(y1,y2,y3)+x3+y3)) ...
    );
G33=@(y1,y2,y3) 1/(16*pi*(1-nu))*( ...
    (3-4*nu)./r1(y1,y2,y3)+(5-12*nu+8*nu^2)./r2(y1,y2,y3)+(x3-y3).^2./r1(y1,y2,y3).^3 ...
    +6*x3.*y3.*(x3+y3).^2./r2(y1,y2,y3).^5+((3-4*nu)*(x3+y3).^2-2*x3.*y3)./r2(y1,y2,y3).^3 ...
    );

% parameterized surface integral, mapping uv space to three-dimensional space
y=@(u,v,A,B,C) A*(1-u)*(1-v)/4+B*(1+u)*(1-v)/4+C*(1+v)/2;

% moment density
m11=lambda*ekk+2*e11;
m12=2*e12;
m13=2*e13;
m22=lambda*ekk+2*e22;
m23=2*e23;
m33=lambda*ekk+2*e33;

% numerical solution with Gauss-Legendre quadrature
legendreGauss=cell(61,1);

% http://keisan.casio.com/exec/system/1329114617

legendreGauss{7}=[...
    -0.94910791234275852452618968404785126240077093767062,	0.12948496616886969327061143267908201832858740225995;
    -0.74153118559939443986386477328078840707414764714139,	0.2797053914892766679014677714237795824869250652266;
    -0.4058451513773971669066064120769614633473820140994,	0.3818300505051189449503697754889751338783650835339;
    0,                                                      0.4179591836734693877551020408163265306122448979592;
    0.40584515137739716690660641207696146334738201409937,	0.38183005050511894495036977548897513387836508353386;
    0.74153118559939443986386477328078840707414764714139,	0.2797053914892766679014677714237795824869250652266;
    0.94910791234275852452618968404785126240077093767062,	0.12948496616886969327061143267908201832858740225995];

legendreGauss{11}=[...
    -0.97822865814605699280393800112285739077142240891978,	0.055668567116173666482753720442548578728515625696898;
    -0.88706259976809529907515776930392726663167575122531,	0.1255803694649046246346942992239401001976157913954;
    -0.73015200557404932409341625203115345804964306202613,	0.1862902109277342514260976414316558916912847480402;
    -0.51909612920681181592572566945860955448022711511993,	0.23319376459199047991852370484317513943179817231696;
    -0.26954315595234497233153198540086152467962186243905,	0.26280454451024666218068886989050919537276467760315;
    0,                                                      0.27292508677790063071448352833634218915604196989478;
    0.26954315595234497233153198540086152467962186243905,	0.26280454451024666218068886989050919537276467760315;
    0.51909612920681181592572566945860955448022711511993,	0.23319376459199047991852370484317513943179817231696;
    0.73015200557404932409341625203115345804964306202613,	0.1862902109277342514260976414316558916912847480402;
    0.88706259976809529907515776930392726663167575122531,	0.1255803694649046246346942992239401001976157913954;
    0.97822865814605699280393800112285739077142240891978,	0.0556685671161736664827537204425485787285156256969];

legendreGauss{15}=[...
    -0.98799251802048542848956571858661258114697281712376,	0.03075324199611726835462839357720441772174814483343;
    -0.9372733924007059043077589477102094712439962735153,	0.07036604748810812470926741645066733846670803275433;
    -0.84820658341042721620064832077421685136625617473699,	0.1071592204671719350118695466858693034155437157581;
    -0.7244177313601700474161860546139380096308992945841,	0.13957067792615431444780479451102832252085027531551;
    -0.57097217260853884753722673725391064123838639628275,	0.16626920581699393355320086048120881113090018009841;
    -0.3941513470775633698972073709810454683627527761587,	0.1861610000155622110268005618664228245062260122779;
    -0.20119409399743452230062830339459620781283645446264,	0.19843148532711157645611832644383932481869255995754;
    0,                                                      0.20257824192556127288062019996751931483866215800948;
    0.20119409399743452230062830339459620781283645446264,	0.19843148532711157645611832644383932481869255995754;
    0.3941513470775633698972073709810454683627527761587,    0.1861610000155622110268005618664228245062260122779;
    0.57097217260853884753722673725391064123838639628275,	0.16626920581699393355320086048120881113090018009841;
    0.7244177313601700474161860546139380096308992945841,    0.13957067792615431444780479451102832252085027531551;
    0.84820658341042721620064832077421685136625617473699,	0.1071592204671719350118695466858693034155437157581;
    0.93727339240070590430775894771020947124399627351531,	0.070366047488108124709267416450667338466708032754331;
    0.98799251802048542848956571858661258114697281712376,	0.03075324199611726835462839357720441772174814483343];

legendreGauss{19}=[ ...
    -0.99240684384358440318901767025326049358931640140321,	0.01946178822972647703631204146443843575290660906929;
    -0.96020815213483003085277884068765152661509150327414,	0.04481422676569960033283815740199421195175422746786;
    -0.90315590361481790164266092853231248780939393405736,	0.06904454273764122658070825800601304496184803168761;
    -0.82271465653714282497892248671271390177453848620683,	0.0914900216224499994644620941238396526609116512966;
    -0.7209661773352293786170958608237816296571418329087,	0.11156664554733399471602390168176599748133185383989;
    -0.60054530466168102346963816494623927986832208273229,	0.12875396253933622767551578485687711705583957709346;
    -0.46457074137596094571726714810410236797628571462414,	0.1426067021736066117757461094419029724756683448245;
    -0.31656409996362983199011732884984491789228521913289,	0.15276604206585966677885540089766299846100826723643;
    -0.16035864564022537586809611574074354950487350047088,	0.15896884339395434764995643946504720167878015819513;
     0,                                                 	0.1610544498487836959791636253209167350399025585785;
     0.16035864564022537586809611574074354950487350047088,	0.1589688433939543476499564394650472016787801581951;
     0.31656409996362983199011732884984491789228521913289,	0.1527660420658596667788554008976629984610082672364;
     0.46457074137596094571726714810410236797628571462414,	0.1426067021736066117757461094419029724756683448245;
     0.60054530466168102346963816494623927986832208273229,	0.1287539625393362276755157848568771170558395770935;
     0.72096617733522937861709586082378162965714183290867,	0.1115666455473339947160239016817659974813318538399;
     0.82271465653714282497892248671271390177453848620683,	0.0914900216224499994644620941238396526609116512966;
     0.90315590361481790164266092853231248780939393405736,	0.06904454273764122658070825800601304496184803168761;
     0.96020815213483003085277884068765152661509150327414,	0.04481422676569960033283815740199421195175422746786;
     0.99240684384358440318901767025326049358931640140321,	0.01946178822972647703631204146443843575290660906929];


legendreGauss{29}=[ ...
    -0.99667944226059658616319153254935388565177345325088,	0.00851690387874640965426381330224980300239888979852;
    -0.98254550526141317487092601578637695610678194890643,	0.01973208505612270598385980164039563114960656819036;
    -0.9572855957780877257982080369808235637375595539541,	0.03074049220209362264440852537461674974711626096871;
    -0.92118023295305878509375343608310642540883934331229,	0.04140206251868283610483001011407692153349078038857;
    -0.87463780492010279041779342125657854691061686311904,	0.0515948269024979239125943811795425979196221106461;
    -0.81818548761525244498957221457878497563002156188229,	0.06120309065707913854210984802390704492406453977913;
    -0.75246285173447713391261007721213819021259196822899,	0.07011793325505127856958148694887917310239488867477;
    -0.67821453760268651515618500539198592638826293001465,	0.07823832713576378382814488865968033136687907695967;
    -0.59628179713822782037958621118898978007818681038558,	0.08547225736617252754534484929720807138169825661359;
    -0.50759295512422764210262791962752015335122384038458,	0.09173775713925876334796641107711080689822902098195;
    -0.41315288817400866389070658603161682332343851873631,	0.0969638340944086063019000748826887591763562779218;
    -0.31403163786763993494819592319104744825844881909083,	0.10109127375991496612182054690749736364756727885253;
    -0.2113522861660010745063757289029374990566508177121,	0.10407331007772937391332847128512006891065222067547;
    -0.10627823013267923017098239243037698091719543929561,	0.1058761550973209414065913278521878930748313137423;
    0,                                                      0.106479381718314244246511126909677568330185031613;
    0.10627823013267923017098239243037698091719543929561,	0.1058761550973209414065913278521878930748313137423;
    0.2113522861660010745063757289029374990566508177121,	0.10407331007772937391332847128512006891065222067547;
    0.31403163786763993494819592319104744825844881909083,	0.10109127375991496612182054690749736364756727885253;
    0.41315288817400866389070658603161682332343851873631,	0.0969638340944086063019000748826887591763562779218;
    0.50759295512422764210262791962752015335122384038458,	0.09173775713925876334796641107711080689822902098195;
    0.5962817971382278203795862111889897800781868103856,	0.08547225736617252754534484929720807138169825661359;
    0.67821453760268651515618500539198592638826293001465,	0.0782383271357637838281448886596803313668790769597;
    0.75246285173447713391261007721213819021259196822899,	0.07011793325505127856958148694887917310239488867477;
    0.81818548761525244498957221457878497563002156188229,	0.06120309065707913854210984802390704492406453977913;
    0.87463780492010279041779342125657854691061686311904,	0.05159482690249792391259438117954259791962211064614;
    0.92118023295305878509375343608310642540883934331229,	0.0414020625186828361048300101140769215334907803886;
    0.9572855957780877257982080369808235637375595539541,	0.03074049220209362264440852537461674974711626096871;
    0.98254550526141317487092601578637695610678194890643,	0.01973208505612270598385980164039563114960656819036;
    0.99667944226059658616319153254935388565177345325088,	0.00851690387874640965426381330224980300239888979852];

legendreGauss{39}=[ ...
    -0.99814738306643290600547230285182044183200974004966,	0.00475294469163510137077621315490694198569508454184;
    -0.99025153685468598363977511724707820470969382104019,	0.01103478893916459424267680545217728945350713937312;
    -0.9760987093334710538448503198895184942907474126801,	0.01725622909372491904080547118335504721371590189046;
    -0.95577521232465227711089189719108220300177179063888,	0.02336938483217816459471234444292469889119755546376;
    -0.9294091484867382296978169643577345618441028095077,	0.0293349559839033785921559863562514256114114671016;
    -0.89716711929299288784829109086081365106569207797562,	0.03511511149813133076106518529723281636666148761973;
    -0.85925293799990615391379743912604237053508444531982,	0.0406732768479338439390565560822614036875724359961;
    -0.81590629743014310435323267840962982772781406052092,	0.045974301108916631884176639399232780598861693034;
    -0.76740124293106349983227240422834294547283411345901,	0.05098466529212940521402103367658133501761308388128;
    -0.7140444358945346791338670361518307896643452818937,	0.05567269034091629990739113978968148270696094807358;
    -0.65617321343201091073442593497625618488010761157697,	0.060008736088596149574941773548819499736678227091846;
    -0.5941534549572779886928900746191060681655915446781,	0.06396538813868238898670640441006364485577691689;
    -0.52837726866043747389634363580822250534071765348543,	0.06751763096623126536302132804644859377515270808348;
    -0.45926051230913604866324663310957236293874915369602,	0.07064300597060876077011493152813113646742968415903;
    -0.38724016397156145585388196563696269601257617271683,	0.07332175341426861738115393286504830644789342425098;
    -0.31277155924818592253599691178562873960190560997754,	0.075536937322836057704784446990809399094187978590022;
    -0.23632551246183576733600632733265044806571922479406,	0.07727455254468201672851163673311729912049525148879;
    -0.15838533999783779992270106136139938446682623686514,	0.0785236132873711767250633009855234715225814764092;
    -0.07944380460875547758191708319264064688500316510674,	0.079276222568368471010155771754507930523804148678266;
    0,                                                   	0.07952762213944285241741819660585099384560677476586;
    0.07944380460875547758191708319264064688500316510674,	0.07927622256836847101015577175450793052380414867827;
    0.1583853399978377999227010613613993844668262368651,	0.0785236132873711767250633009855234715225814764092;
    0.23632551246183576733600632733265044806571922479406,	0.0772745525446820167285116367331172991204952514888;
    0.3127715592481859225359969117856287396019056099775,	0.07553693732283605770478444699080939909418797859;
    0.38724016397156145585388196563696269601257617271683,	0.073321753414268617381153932865048306447893424251;
    0.45926051230913604866324663310957236293874915369602,	0.070643005970608760770114931528131136467429684159028;
    0.52837726866043747389634363580822250534071765348543,	0.0675176309662312653630213280464485937751527080835;
    0.59415345495727798869289007461910606816559154467813,	0.06396538813868238898670640441006364485577691688998;
    0.65617321343201091073442593497625618488010761157697,	0.06000873608859614957494177354881949973667822709185;
    0.7140444358945346791338670361518307896643452818937,	0.05567269034091629990739113978968148270696094807358;
    0.76740124293106349983227240422834294547283411345901,	0.0509846652921294052140210336765813350176130838813;
    0.81590629743014310435323267840962982772781406052092,	0.045974301108916631884176639399232780598861693034;
    0.85925293799990615391379743912604237053508444531982,	0.0406732768479338439390565560822614036875724359961;
    0.89716711929299288784829109086081365106569207797562,	0.03511511149813133076106518529723281636666148761973;
    0.9294091484867382296978169643577345618441028095077,	0.0293349559839033785921559863562514256114114671016;
    0.95577521232465227711089189719108220300177179063888,	0.0233693848321781645947123444429246988911975554638;
    0.9760987093334710538448503198895184942907474126801,	0.01725622909372491904080547118335504721371590189046;
    0.99025153685468598363977511724707820470969382104019,	0.01103478893916459424267680545217728945350713937312;
    0.99814738306643290600547230285182044183200974004966,	0.00475294469163510137077621315490694198569508454184];


legendreGauss{48}=[ ...
	-0.99877100725242611860054149156311364008893765027672,	0.00315334605230583863267731154389148757828393883169;
    -0.99353017226635075754792875084907411835661474959467,	0.0073275539012762621023839796217865500587079025592;
    -0.9841245837228268577445836000265988305892392234174,	0.01147723457923453948959266760909162808642050630875;
    -0.9705915925462472504614119838006600573024339116309,	0.01557931572294384872817695583446031397637626899155;
    -0.95298770316043086072296066602571834320854133182392,	0.0196161604573555278144607196522127096958130377341;
    -0.93138669070655433311417438010160126771999708561895,	0.02357076083932437914051930137844923022172973852219;
    -0.90587913671556967282207483567101178831226219982741,	0.02742650970835694820007383626250582045118415516165;
    -0.87657202027424788590569355480509675456164853372996,	0.03116722783279808890206575684635441945428534148357;
    -0.84358826162439353071108984451965604987088701173755,	0.034777222564770438892548585963802410597281396907068;
    -0.80706620402944262708255304302453844597301302946042,	0.038241351065830706317217256523715617863823968355;
    -0.76715903251574033925385543752296905362264233084821,	0.041545082943464749214058822361064797753472826034038;
    -0.72403413092381465467448223349366524658509281228072,	0.0446745608566942804194485871258503949884627868625;
    -0.677872379632663905211851280675909058849954679026,	0.04761665849249047482590662347892983015799806674345;
    -0.62886739677651362399516493306999465202490899979016,	0.05035903555385447495780761908786560603299409302591;
    -0.5772247260839727038178092385404787728539972861402,	0.052890189485193667095505056264698914661726485633109;
    -0.52316097472223303367822586913750852628918762181188,	0.0551995036999841628682034951916354390044509256076;
    -0.46690290475095840454492886165079850923681210425852,	0.05727729210040321570515023468470057624152712300411;
    -0.40868648199071672991622549581463328645992284299489,	0.059114839698395635746474817433519910659655602557055;
    -0.3487558862921607381598179372704079161343096499684,	0.06070443916589388005296923202782047788526086425648;
    -0.28736248735545557673588646131679768785155830580104,	0.062039423159892663904197784137598518306383399665092;
    -0.22476379039468906122486544017469227743856180404166,	0.0631141922862540256571260227502333181274136433711;
    -0.1612223560688917180564373907834976947743743797419,	0.063924238584648186623906201825515408918974084982643;
    -0.09700469920946269893005395585362452015273622930094,	0.0644661644359500822065041936577050657256919244555;
    -0.03238017096286936203332224315213444204596280236152,	0.06473769681268392250302493873659155355208191894664;
     0.03238017096286936203332224315213444204596280236152,	0.06473769681268392250302493873659155355208191894664;
     0.09700469920946269893005395585362452015273622930094,	0.0644661644359500822065041936577050657256919244555;
     0.1612223560688917180564373907834976947743743797419,	0.0639242385846481866239062018255154089189740849826;
     0.22476379039468906122486544017469227743856180404166,	0.0631141922862540256571260227502333181274136433711;
     0.28736248735545557673588646131679768785155830580104,	0.06203942315989266390419778413759851830638339966509;
     0.34875588629216073815981793727040791613430964996839,	0.0607044391658938800529692320278204778852608642565;
     0.40868648199071672991622549581463328645992284299489,	0.05911483969839563574647481743351991065965560255706;
     0.46690290475095840454492886165079850923681210425852,	0.0572772921004032157051502346847005762415271230041;
     0.52316097472223303367822586913750852628918762181188,	0.0551995036999841628682034951916354390044509256076;
     0.5772247260839727038178092385404787728539972861402,	0.0528901894851936670955050562646989146617264856331;
     0.62886739677651362399516493306999465202490899979016,	0.05035903555385447495780761908786560603299409302591;
     0.67787237963266390521185128067590905884995467902605,	0.04761665849249047482590662347892983015799806674345;
     0.72403413092381465467448223349366524658509281228072,	0.0446745608566942804194485871258503949884627868625;
     0.76715903251574033925385543752296905362264233084821,	0.041545082943464749214058822361064797753472826034;
     0.80706620402944262708255304302453844597301302946042,	0.038241351065830706317217256523715617863823968355;
     0.8435882616243935307110898445196560498708870117376,	0.03477722256477043889254858596380241059728139690707;
     0.87657202027424788590569355480509675456164853373, 	0.0311672278327980889020657568463544194542853414836;
     0.90587913671556967282207483567101178831226219982741,	0.0274265097083569482000738362625058204511841551617;
     0.93138669070655433311417438010160126771999708561895,	0.02357076083932437914051930137844923022172973852219;
     0.95298770316043086072296066602571834320854133182392,	0.01961616045735552781446071965221270969581303773413;
     0.97059159254624725046141198380066005730243391163088,	0.01557931572294384872817695583446031397637626899155;
     0.9841245837228268577445836000265988305892392234174,	0.01147723457923453948959266760909162808642050630875;
     0.99353017226635075754792875084907411835661474959467,	0.0073275539012762621023839796217865500587079025592;
     0.9987710072524261186005414915631136400889376502767,	0.00315334605230583863267731154389148757828393883169];

legendreGauss{61}=[ ...
    -0.999235597631363471731862259569130862853442048503,	0.00196145336167028267177306446192684809414929406187;
    -0.99597459981512023426801276071228827628684886478913,	0.00456092400601241718453893257974278827081304337832;
    -0.9901167452325170509655316856999165016287567811765,	0.00715235499174908958583399309169417557127328166703;
    -0.98167601128403707968517254075131942216007406209493,	0.00972546183035613373613667175557011832684334637992;
    -0.97067425883318290824740845354791027256833637675753,	0.01227326350781210462927780035952613135449266049143;
    -0.95714015191298409137208028742121526110444034092702,	0.01478906588493791454617851578148826126037152853378;
    -0.94110898668136114747754386097778900828047471105846,	0.01726629298761374359443321277243662118418805452052;
    -0.92262258138295526125755150055584544545063516129872,	0.01969847774610118133051783021268456666579663613604;
    -0.90172916247400117064204004268351116465311752108726,	0.02207927314831904400247576562818162743617212408293;
    -0.87848323721488103247894751389335263145955452028042,	0.02440246718754420291534065973246753232577019850402;
    -0.85294545084766344556484690467417053833995325998446,	0.02666199852415088966281070684496801956405541912412;
    -0.82518242810865995066428189245084752700639465995958,	0.028851972088183401504341775545590252543507119934694;
    -0.79526659928235964915204880275254626487323687020583,	0.03096667436839739482469794661248938376670927102753;
    -0.76327601117231219714591466849435037104451915109862,	0.03300058827590741063272368344100449280740893980815;
    -0.72929412344946510968895669013553179857637101453818,	0.03494840751653335109085191920921927825500947094151;
    -0.69340959089449115549918436396023173728437756854619,	0.03680505042315481738432104286406746087590072359162;
    -0.65571603209507087169918577783630688700084734644168,	0.038565673207008172746152047629653802829222933215;
    -0.61631178519792172470961685932859103589567168396123,	0.0402256825909982473676399837575115791559285285247;
    -0.57529965135083061860036998194497287899733192533773,	0.041780747790888492066675572351073007512697374083565;
    -0.53278662650292526563848171637307164032481991094547,	0.0432268118124960979010436458226828733502967443538;
    -0.488883622262252118820698511425829597516530924725,	0.04456010203508348827154141983108835293335633519477;
    -0.4437051765385316019955893010658225593823246374722,	0.04577714005314595937133983384682479196153481455209;
    -0.3973691547257566091782918473685892218128728926937,	0.0468747507508090659764294455885907359337142723544;
    -0.3499964422040668345334344770620938497875046174049,	0.04785007058509560716183342751869966271340125297024;
    -0.30171062896303071260448652545697712524119763287, 	0.04870055505641152608753008831470560675869837218145;
    -0.25263768716905349583369086333440510622751074756006,	0.049423985346735589939968776651116499416659656667749;
    -0.20290564251805849922694720334305705450143905185836,	0.05001847410817825342505161500635862544205591923667;
    -0.15264424023081530052950676177347741688353897661746,	0.05048247038679740464814446518815126917289915200009;
    -0.10198460656227406895720840476436564197678981705575,	0.05081476366881834320770052922347869734784234385657;
    -0.05105890670797434936688750061890079807018425404126,	0.0510144870386972635437350805738520107682023754234;
    0,                                                      0.0510811194407862179779210956063098528020576264041;
    0.05105890670797434936688750061890079807018425404126,	0.0510144870386972635437350805738520107682023754234;
    0.1019846065622740689572084047643656419767898170558,	0.05081476366881834320770052922347869734784234385657;
    0.15264424023081530052950676177347741688353897661746,	0.05048247038679740464814446518815126917289915200009;
    0.2029056425180584992269472033430570545014390518584,	0.0500184741081782534250516150063586254420559192367;
    0.25263768716905349583369086333440510622751074756006,	0.04942398534673558993996877665111649941665965666775;
    0.30171062896303071260448652545697712524119763287,  	0.04870055505641152608753008831470560675869837218145;
    0.34999644220406683453343447706209384978750461740488,	0.0478500705850956071618334275186996627134012529702;
    0.39736915472575660917829184736858922181287289269372,	0.0468747507508090659764294455885907359337142723544;
    0.44370517653853160199558930106582255938232463747224,	0.04577714005314595937133983384682479196153481455209;
    0.488883622262252118820698511425829597516530924725, 	0.044560102035083488271541419831088352933356335194765;
    0.53278662650292526563848171637307164032481991094547,	0.0432268118124960979010436458226828733502967443538;
    0.57529965135083061860036998194497287899733192533773,	0.04178074779088849206667557235107300751269737408356;
    0.61631178519792172470961685932859103589567168396123,	0.04022568259099824736763998375751157915592852852474;
    0.65571603209507087169918577783630688700084734644168,	0.03856567320700817274615204762965380282922293321496;
    0.69340959089449115549918436396023173728437756854619,	0.03680505042315481738432104286406746087590072359162;
    0.72929412344946510968895669013553179857637101453818,	0.0349484075165333510908519192092192782550094709415;
    0.76327601117231219714591466849435037104451915109862,	0.0330005882759074106327236834410044928074089398082;
    0.79526659928235964915204880275254626487323687020583,	0.03096667436839739482469794661248938376670927102753;
    0.82518242810865995066428189245084752700639465995958,	0.0288519720881834015043417755455902525435071199347;
    0.85294545084766344556484690467417053833995325998446,	0.026661998524150889662810706844968019564055419124121;
    0.87848323721488103247894751389335263145955452028042,	0.024402467187544202915340659732467532325770198504;
    0.90172916247400117064204004268351116465311752108726,	0.02207927314831904400247576562818162743617212408293;
    0.92262258138295526125755150055584544545063516129872,	0.01969847774610118133051783021268456666579663613604;
    0.94110898668136114747754386097778900828047471105846,	0.01726629298761374359443321277243662118418805452052;
    0.95714015191298409137208028742121526110444034092702,	0.0147890658849379145461785157814882612603715285338;
    0.97067425883318290824740845354791027256833637675753,	0.01227326350781210462927780035952613135449266049143;
    0.98167601128403707968517254075131942216007406209493,	0.0097254618303561337361366717555701183268433463799;
    0.9901167452325170509655316856999165016287567811765,	0.00715235499174908958583399309169417557127328166703;
    0.99597459981512023426801276071228827628684886478913,	0.00456092400601241718453893257974278827081304337832;
    0.999235597631363471731862259569130862853442048503, 	0.00196145336167028267177306446192684809414929406187];
 
x=legendreGauss{optionStruct.n}(:,1);
w=legendreGauss{optionStruct.n}(:,2);

u1=zeros(size(x1));
u2=zeros(size(x2));
u3=zeros(size(x3));

for k=1:length(x)
    for j=1:length(x)
        u1=u1+w(j)*w(k)*(1-x(k))*IU1(x(j),x(k));
        u2=u2+w(j)*w(k)*(1-x(k))*IU2(x(j),x(k));
        u3=u3+w(j)*w(k)*(1-x(k))*IU3(x(j),x(k));
    end
end

    function d = IU1(u,v)
        % function IU1 is the integrand for displacement component u1
        d=zeros(size(x1));
        if (m11*nD(1)+m12*nD(2)+m13*nD(3)) ~= 0
            d=d+ABC/4*(m11*nD(1)+m12*nD(2)+m13*nD(3))*G11(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m12*nD(1)+m22*nD(2)+m23*nD(3)) ~= 0
            d=d+ABC/4*(m12*nD(1)+m22*nD(2)+m23*nD(3))*G21(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m13*nD(1)+m23*nD(2)+m33*nD(3)) ~= 0
            d=d+ABC/4*(m13*nD(1)+m23*nD(2)+m33*nD(3))*G31(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m11*nA(1)+m12*nA(2)+m13*nA(3)) ~= 0
            d=d+BCD/4*(m11*nA(1)+m12*nA(2)+m13*nA(3))*G11(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m12*nA(1)+m22*nA(2)+m23*nA(3)) ~= 0
            d=d+BCD/4*(m12*nA(1)+m22*nA(2)+m23*nA(3))*G21(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m13*nA(1)+m23*nA(2)+m33*nA(3)) ~= 0
            d=d+BCD/4*(m13*nA(1)+m23*nA(2)+m33*nA(3))*G31(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m11*nB(1)+m12*nB(2)+m13*nB(3)) ~= 0
            d=d+CDA/4*(m11*nB(1)+m12*nB(2)+m13*nB(3))*G11(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m12*nB(1)+m22*nB(2)+m23*nB(3)) ~= 0
            d=d+CDA/4*(m12*nB(1)+m22*nB(2)+m23*nB(3))*G21(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m13*nB(1)+m23*nB(2)+m33*nB(3)) ~= 0
            d=d+CDA/4*(m13*nB(1)+m23*nB(2)+m33*nB(3))*G31(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m11*nC(1)+m12*nC(2)+m13*nC(3)) ~= 0
            d=d+DAB/4*(m11*nC(1)+m12*nC(2)+m13*nC(3))*G11(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
        if (m12*nC(1)+m22*nC(2)+m23*nC(3)) ~= 0
            d=d+DAB/4*(m12*nC(1)+m22*nC(2)+m23*nC(3))*G21(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
        if (m13*nC(1)+m23*nC(2)+m33*nC(3)) ~= 0
            d=d+DAB/4*(m13*nC(1)+m23*nC(2)+m33*nC(3))*G31(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
    end

    function d = IU2(u,v)
        % function IU2 is the integrand for displacement component u2
        d=zeros(size(x2));
        if (m11*nD(1)+m12*nD(2)+m13*nD(3)) ~= 0
            d=d+ABC/4*(m11*nD(1)+m12*nD(2)+m13*nD(3))*G12(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m12*nD(1)+m22*nD(2)+m23*nD(3)) ~= 0
            d=d+ABC/4*(m12*nD(1)+m22*nD(2)+m23*nD(3))*G22(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m13*nD(1)+m23*nD(2)+m33*nD(3)) ~= 0
            d=d+ABC/4*(m13*nD(1)+m23*nD(2)+m33*nD(3))*G32(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m11*nA(1)+m12*nA(2)+m13*nA(3)) ~= 0
            d=d+BCD/4*(m11*nA(1)+m12*nA(2)+m13*nA(3))*G12(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m12*nA(1)+m22*nA(2)+m23*nA(3)) ~= 0
            d=d+BCD/4*(m12*nA(1)+m22*nA(2)+m23*nA(3))*G22(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m13*nA(1)+m23*nA(2)+m33*nA(3)) ~= 0
            d=d+BCD/4*(m13*nA(1)+m23*nA(2)+m33*nA(3))*G32(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m11*nB(1)+m12*nB(2)+m13*nB(3)) ~= 0
            d=d+CDA/4*(m11*nB(1)+m12*nB(2)+m13*nB(3))*G12(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m12*nB(1)+m22*nB(2)+m23*nB(3)) ~= 0
            d=d+CDA/4*(m12*nB(1)+m22*nB(2)+m23*nB(3))*G22(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m13*nB(1)+m23*nB(2)+m33*nB(3)) ~= 0
            d=d+CDA/4*(m13*nB(1)+m23*nB(2)+m33*nB(3))*G32(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m11*nC(1)+m12*nC(2)+m13*nC(3)) ~= 0
            d=d+DAB/4*(m11*nC(1)+m12*nC(2)+m13*nC(3))*G12(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
        if (m12*nC(1)+m22*nC(2)+m23*nC(3)) ~= 0
            d=d+DAB/4*(m12*nC(1)+m22*nC(2)+m23*nC(3))*G22(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
        if (m13*nC(1)+m23*nC(2)+m33*nC(3)) ~= 0
            d=d+DAB/4*(m13*nC(1)+m23*nC(2)+m33*nC(3))*G32(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
    end

    function d = IU3(u,v)
        % function IU3 is the integrand for displacement component u3
        d=zeros(size(x3));
        if (m11*nD(1)+m12*nD(2)+m13*nD(3)) ~= 0
            d=d+ABC/4*(m11*nD(1)+m12*nD(2)+m13*nD(3))*G13(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m12*nD(1)+m22*nD(2)+m23*nD(3)) ~= 0
            d=d+ABC/4*(m12*nD(1)+m22*nD(2)+m23*nD(3))*G23(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m13*nD(1)+m23*nD(2)+m33*nD(3)) ~= 0
            d=d+ABC/4*(m13*nD(1)+m23*nD(2)+m33*nD(3))*G33(y(u,v,A(1),B(1),C(1)),y(u,v,A(2),B(2),C(2)),y(u,v,A(3),B(3),C(3)));
        end
        if (m11*nA(1)+m12*nA(2)+m13*nA(3)) ~= 0
            d=d+BCD/4*(m11*nA(1)+m12*nA(2)+m13*nA(3))*G13(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m12*nA(1)+m22*nA(2)+m23*nA(3)) ~= 0
            d=d+BCD/4*(m12*nA(1)+m22*nA(2)+m23*nA(3))*G23(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m13*nA(1)+m23*nA(2)+m33*nA(3)) ~= 0
            d=d+BCD/4*(m13*nA(1)+m23*nA(2)+m33*nA(3))*G33(y(u,v,B(1),C(1),D(1)),y(u,v,B(2),C(2),D(2)),y(u,v,B(3),C(3),D(3)));
        end
        if (m11*nB(1)+m12*nB(2)+m13*nB(3)) ~= 0
            d=d+CDA/4*(m11*nB(1)+m12*nB(2)+m13*nB(3))*G13(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m12*nB(1)+m22*nB(2)+m23*nB(3)) ~= 0
            d=d+CDA/4*(m12*nB(1)+m22*nB(2)+m23*nB(3))*G23(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m13*nB(1)+m23*nB(2)+m33*nB(3)) ~= 0
            d=d+CDA/4*(m13*nB(1)+m23*nB(2)+m33*nB(3))*G33(y(u,v,C(1),D(1),A(1)),y(u,v,C(2),D(2),A(2)),y(u,v,C(3),D(3),A(3)));
        end
        if (m11*nC(1)+m12*nC(2)+m13*nC(3)) ~= 0
            d=d+DAB/4*(m11*nC(1)+m12*nC(2)+m13*nC(3))*G13(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
        if (m12*nC(1)+m22*nC(2)+m23*nC(3)) ~= 0
            d=d+DAB/4*(m12*nC(1)+m22*nC(2)+m23*nC(3))*G23(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
        if (m13*nC(1)+m23*nC(2)+m33*nC(3)) ~= 0
            d=d+DAB/4*(m13*nC(1)+m23*nC(2)+m33*nC(3))*G33(y(u,v,D(1),A(1),B(1)),y(u,v,D(2),A(2),B(2)),y(u,v,D(3),A(3),B(3)));
        end
    end

end

function p = validateN(x)
if ~(x > 0)
    error('MATLAB:invalid','invalid number of integration points');
end
p = true;
end

