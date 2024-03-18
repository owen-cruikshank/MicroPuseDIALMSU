%model of returns

clear all
%Pulse length model
h = 6.62607004e-34;
c = 2.9989e8;   % Speed of light m/s

%Input Variables
dt = .5;         %pulse duration in microseconds
nt = 40;        %number of points to calculate in each dt. THIS MUST BE AN EVEN NUMBER
nt = 2;
nt = 20;
nt = 85;

Ts = 273;       %Surface Temperature in K      
Ps = .825;         %Surface Pressure in atm

Ts = 288.15;
Ps = 1.01325e3*0.000986923;

%lambda_online = 769.7958;
%nu_online = 10^7./lambda_online;                    %[cm-1] Online wavenumber
WV=0;

% Range Vector
dr = c*dt*1e-6/nt;
rm = 0.1:dr:18000;             %Range vector in m
rkm =rm./1000;               %Range vector in km
size_r = length(rm);

% Temperature and Pressure Profiles
T = Ts-rkm.*6.5;           %Temperature profile in K
P = Ps.*(Ts./T).^-5.2199;   %Pressure profile in atm
T = T';
P = P';

%Radiosonde proflies
sondepath = 'C:\Users\Owen\OneDrive - Montana State University\Research\O2 DIAL\Data\MSU data\Radiosondes';
span_days = datetime(2022,6,22,'TimeZone','UTC');%yyyy,mm,dd)
[sonde_datetime,sondeStruc] =  COBradiosonde(sondepath,span_days);

sondeT = sondeStruc(1).T;
sondeP = sondeStruc(1).P.*0.000986923;
sondeH = sondeStruc(1).Height-sondeStruc(1).Height(1);

% % T = interp1(sondeH,sondeT,rm)';
% % P = interp1(sondeH,sondeP,rm)';

A = pi*.406.^2-pi*.2.^2;

eta_D = .6;%detector
eta_R = .0001;

E_pulse_on = 43 * 10^-6;%pulse energy (J)
E_pulse_off = 43 * 10^-6;%pulse energy (J)

%--Outgoint
pulse_rate = 7000;%(Hz)
avg_time = 10*60*60;%(sec)

%%

%---Backscatter
%%%Bm = 374280*101.325*P./T./(lambda_online.^4);
% 
% Ba = Bm;
% 
% Ba = zeros(size(Bm)).*0.01;
% Ba(1:512) = 3.*Bm(1:512) ;
% %Ba(512:600) = 2.*Bm(512:600) -2.*Bm(512:600).*(1:89)'.*(1/89);
% 
% Ba(512:700) = 3.*Bm(512:700) -3.*Bm(512:700).*(1:189)'.*(1/189);
%%
BSR=[
NaN
NaN
NaN
NaN
NaN
3.33348686081786
3.36551751422900
3.24830120190958
3.51397407672946
3.24886576227332
3.28256129485638
3.20222384944337
2.99930868312942
2.86916198265844
2.64385031814741
2.47578061846368
2.42595882010490
2.36780471865549
2.44772042597672
2.37986315656486
2.40141244360314
2.52147184604419
2.65145044377840
2.63842063670806
2.70310582286966
2.72762629409563
2.67901534569066
2.55623055902631
2.53993907663394
2.49624738449374
2.69396749176289
2.90236719722746
3.13748868430628
3.10623558305595
3.28708912417195
3.02949145827879
2.40922813682599
1.76943604070225
1.41563749633709
1.19546276867337
1.10688539707404
1.11987369923108
1.09801303504226
1.13843123587360
1.11719262951284
1.10561191716656
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN];

Ba=[
    NaN
NaN
NaN
NaN
NaN
6.32907992830923e-07
6.37710590617329e-07
6.02431327959596e-07
6.69518780572237e-07
5.95261370654154e-07
6.00484865622041e-07
5.75797461542331e-07
5.19529228403471e-07
4.82716417259260e-07
4.21905813186831e-07
3.76423065828719e-07
3.61456164568666e-07
3.44556130881516e-07
3.62410512957459e-07
3.43261612072238e-07
3.46434565652825e-07
3.73747293067504e-07
4.03117256846504e-07
3.97407212971547e-07
4.10477358414255e-07
4.13739821262410e-07
3.99534967784686e-07
3.67950426831804e-07
3.61765163599222e-07
3.49242422286174e-07
3.92845292441800e-07
4.38324953152219e-07
4.89309294805676e-07
4.79023446621442e-07
5.16767802649679e-07
4.55569274199370e-07
3.14264704622478e-07
1.70461440408553e-07
9.14744603992016e-08
4.27339292735071e-08
2.32136212586037e-08
2.58616491946922e-08
2.10046738413474e-08
2.94684951760562e-08
2.47803728997268e-08
2.21817605779635e-08
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN
NaN];
Bm=[
2.79528095465145e-07
2.77856427388800e-07
2.76190633539672e-07
2.74530707835691e-07
2.72876644185929e-07
2.71228436490573e-07
2.69586078640885e-07
2.67949564519169e-07
2.66318887998735e-07
2.64694042943861e-07
2.63075023209760e-07
2.61461822642539e-07
2.59854435079168e-07
2.58252854347438e-07
2.56657074265928e-07
2.55067088643964e-07
2.53482891281583e-07
2.51904475969496e-07
2.50331836489049e-07
2.48764966612182e-07
2.47203860101395e-07
2.45648510709708e-07
2.44098912180616e-07
2.42555058248059e-07
2.41016942636376e-07
2.39484559060264e-07
2.37957901224742e-07
2.36436962825110e-07
2.34921737546905e-07
2.33412219065860e-07
2.31908401047868e-07
2.30410277148934e-07
2.28917841015135e-07
2.27431086282582e-07
2.25950006577369e-07
2.24474595515538e-07
2.23004846703032e-07
2.21540753735652e-07
2.20082310199014e-07
2.18629509668504e-07
2.17182345709234e-07
2.15740811875998e-07
2.14304901713227e-07
2.12874608754940e-07
2.11449926524705e-07
2.10030848535588e-07
2.08617368290107e-07
2.07209479280188e-07
2.05807174987115e-07
2.04410448881485e-07
2.03019294423160e-07
2.01633705061216e-07
2.00253674233901e-07
1.98879195368578e-07
1.97510261881683e-07
1.96146867178671e-07
1.94789004653969e-07
1.93436667690925e-07
1.92089849661757e-07
1.90748543927500e-07
1.89412743837961e-07
1.88082442731659e-07
1.86757633935780e-07
1.85438310766120e-07
1.84124466527035e-07
1.82816094511386e-07
1.81513188000485e-07
1.80215740264043e-07
1.78923744560112e-07
1.77637194135034e-07
1.76356082223384e-07
1.75080402047912e-07
1.73810146819491e-07
1.72545309737057e-07
1.71285883987552e-07
1.70031862745871e-07
1.68783239174797e-07
1.67540006424946e-07
1.66302157634712e-07
1.65069685930198e-07];

hsrlRange=[
    22.4221717500000
97.3702862500000
172.318400750000
247.266515250000
322.214629750000
397.162744250000
472.110858750000
547.058973250000
622.007087750000
696.955202250000
771.903316750000
846.851431250000
921.799545750000
996.747660250000
1071.69577475000
1146.64388925000
1221.59200375000
1296.54011825000
1371.48823275000
1446.43634725000
1521.38446175000
1596.33257625000
1671.28069075000
1746.22880525000
1821.17691975000
1896.12503425000
1971.07314875000
2046.02126325000
2120.96937775000
2195.91749225000
2270.86560675000
2345.81372125000
2420.76183575000
2495.70995025000
2570.65806475000
2645.60617925000
2720.55429375000
2795.50240825000
2870.45052275000
2945.39863725000
3020.34675175000
3095.29486625000
3170.24298075000
3245.19109525000
3320.13920975000
3395.08732425000
3470.03543875000
3544.98355325000
3619.93166775000
3694.87978225000
3769.82789675000
3844.77601125000
3919.72412575000
3994.67224025000
4069.62035475000
4144.56846925000
4219.51658375000
4294.46469825000
4369.41281275000
4444.36092725000
4519.30904175000
4594.25715625000
4669.20527075000
4744.15338525000
4819.10149975000
4894.04961425000
4968.99772875000
5043.94584325000
5118.89395775000
5193.84207225000
5268.79018675000
5343.73830125000
5418.68641575000
5493.63453025000
5568.58264475000
5643.53075925000
5718.47887375000
5793.42698825000
5868.37510275000
5943.32321725000];

%%
% Ba=[4.6613204e-07
% 4.1047431e-07
% 5.7582895e-07
% 1.3570170e-07
% 1.3161906e-07
% 1.2765082e-07
% 1.2276352e-07
% 1.2303390e-07
% 1.1718502e-07
% 1.1898334e-07
% 1.0745936e-07
% 1.1723892e-07
% 1.1649375e-07
% 1.1805292e-07
% 1.1192841e-07
% 1.0862890e-07
% 1.1089907e-07
% 1.1256763e-07
% 1.0898042e-07
% 1.1106663e-07
% 1.0765420e-07
% 1.1204284e-07
% 1.0690287e-07
% 1.1209287e-07
% 1.1080180e-07
% 1.1635665e-07
% 1.1331431e-07
% 1.1385733e-07
% 1.1213751e-07
% 1.1782267e-07
% 1.1645337e-07
% 1.1585025e-07
% 1.2495643e-07
% 1.1891133e-07
% 1.1513252e-07
% 1.2463370e-07
% 1.2590860e-07
% 1.3135123e-07
% 1.2629938e-07
% 1.2727517e-07
% 1.2929412e-07
% 1.2085943e-07
% 1.0104966e-07
% 8.8880249e-08
% 6.3357447e-08
% 3.6959669e-08
% 2.3852035e-08
% 1.8021296e-08
% 2.4797666e-08
% 3.0249964e-08
% 2.7557578e-08
% 2.0404345e-08
% 3.3469462e-08
% 1.9044329e-08
% 2.3232468e-08
% 2.5082905e-08
% 3.0843974e-08
% 3.8060968e-08
% 3.2537102e-08
% 2.0391640e-08
% 3.5210391e-08
% 5.2044104e-08
% 3.6961101e-08
% 5.1240381e-08
% 2.5100039e-08
% 2.6902645e-08
% 4.0523894e-08
% 5.5383019e-08
% 4.0988994e-08
% 5.3440385e-08
% 3.7596585e-08
% 3.8991800e-08
% 4.8220500e-08
% 7.6745209e-08
% 3.6774896e-08
% 4.2631818e-08
% 4.2929582e-08
% 4.7819249e-08
% 5.1587488e-08
% 4.0568914e-08
% 4.7318927e-08];
% Bm = [
% 2.81186439741600e-07
% 2.79263655924774e-07
% 2.77350920847412e-07
% 2.75448198244344e-07
% 2.73555451922439e-07
% 2.71672645760581e-07
% 2.69799743709638e-07
% 2.67936709792442e-07
% 2.66083508103758e-07
% 2.64240102810260e-07
% 2.62406458150507e-07
% 2.60582538434914e-07
% 2.58768308045728e-07
% 2.56963731436998e-07
% 2.55168773134555e-07
% 2.53383397735982e-07
% 2.51607569910587e-07
% 2.49841254399380e-07
% 2.48084416015045e-07
% 2.46337019641911e-07
% 2.44599030235931e-07
% 2.42870412824651e-07
% 2.41151132507188e-07
% 2.39441154454198e-07
% 2.37740443907855e-07
% 2.36048966181821e-07
% 2.34366686661220e-07
% 2.32693570802613e-07
% 2.31029584133969e-07
% 2.29374692254641e-07
% 2.27728860835338e-07
% 2.26092055618096e-07
% 2.24464242416257e-07
% 2.22845387114436e-07
% 2.21235455668496e-07
% 2.19634414105526e-07
% 2.18042228523805e-07
% 2.16458865092785e-07
% 2.14884290053056e-07
% 2.13318469716322e-07
% 2.11761370465376e-07
% 2.10212958754071e-07
% 2.08673201107289e-07
% 2.07142064120924e-07
% 2.05619514461842e-07
% 2.04105518867866e-07
% 2.02600044147738e-07
% 2.01103057181101e-07
% 1.99614524918464e-07
% 1.98134414381180e-07
% 1.96662692661416e-07
% 1.95199326922126e-07
% 1.93744284397023e-07
% 1.92297532390554e-07
% 1.90859038277867e-07
% 1.89428769504791e-07
% 1.88006693587801e-07
% 1.86592778113996e-07
% 1.85186990741066e-07
% 1.83789299197269e-07
% 1.82399671281402e-07
% 1.81018074862770e-07
% 1.79644477881162e-07
% 1.78278848346821e-07
% 1.76921154340418e-07
% 1.75571364013022e-07
% 1.74229445586070e-07
% 1.72895367351346e-07
% 1.71569097670945e-07
% 1.70250604977249e-07
% 1.68939857772899e-07
% 1.67636824630766e-07
% 1.66341474193921e-07
% 1.65053775175609e-07
% 1.63773696359219e-07
% 1.62501206598258e-07
% 1.61236274816321e-07
% 1.59978870007061e-07
% 1.58728961234162e-07
% 1.57486517631313e-07
% 1.56251508402174e-07];
% hsrlRange = [
%     0
% 74.9481145000000
% 149.896229000000
% 224.844343500000
% 299.792458000000
% 374.740572500000
% 449.688687000000
% 524.636801500000
% 599.584916000000
% 674.533030500000
% 749.481145000000
% 824.429259500000
% 899.377374000000
% 974.325488500000
% 1049.27360300000
% 1124.22171750000
% 1199.16983200000
% 1274.11794650000
% 1349.06606100000
% 1424.01417550000
% 1498.96229000000
% 1573.91040450000
% 1648.85851900000
% 1723.80663350000
% 1798.75474800000
% 1873.70286250000
% 1948.65097700000
% 2023.59909150000
% 2098.54720600000
% 2173.49532050000
% 2248.44343500000
% 2323.39154950000
% 2398.33966400000
% 2473.28777850000
% 2548.23589300000
% 2623.18400750000
% 2698.13212200000
% 2773.08023650000
% 2848.02835100000
% 2922.97646550000
% 2997.92458000000
% 3072.87269450000
% 3147.82080900000
% 3222.76892350000
% 3297.71703800000
% 3372.66515250000
% 3447.61326700000
% 3522.56138150000
% 3597.50949600000
% 3672.45761050000
% 3747.40572500000
% 3822.35383950000
% 3897.30195400000
% 3972.25006850000
% 4047.19818300000
% 4122.14629750000
% 4197.09441200000
% 4272.04252650000
% 4346.99064100000
% 4421.93875550000
% 4496.88687000000
% 4571.83498450000
% 4646.78309900000
% 4721.73121350000
% 4796.67932800000
% 4871.62744250000
% 4946.57555700000
% 5021.52367150000
% 5096.47178600000
% 5171.41990050000
% 5246.36801500000
% 5321.31612950000
% 5396.26424400000
% 5471.21235850000
% 5546.16047300000
% 5621.10858750000
% 5696.05670200000
% 5771.00481650000
% 5845.95293100000
% 5920.90104550000
% 5995.84916000000];
% %
% Ba = [
%     5.2134612e-07
% 4.3726834e-07
% 3.5601121e-07
% 2.8319079e-07
% 1.7023110e-07
% 1.6213863e-07
% 1.5790022e-07
% 1.5479763e-07
% 1.5162421e-07
% 1.5018401e-07
% 1.4886989e-07
% 1.4750792e-07
% 1.4606844e-07
% 1.4507161e-07
% 1.4369343e-07
% 1.4280425e-07
% 1.4168602e-07
% 1.4049969e-07
% 1.3887345e-07
% 1.3709148e-07
% 1.3566243e-07
% 1.3395567e-07
% 1.3284908e-07
% 1.3128975e-07
% 1.3008751e-07
% 1.2933755e-07
% 1.2791492e-07
% 1.2812750e-07
% 1.2818411e-07
% 1.2834947e-07
% 1.2861798e-07
% 1.2808525e-07
% 1.2714708e-07
% 1.2632697e-07
% 1.2509165e-07
% 1.2426499e-07
% 1.2363002e-07
% 1.2180534e-07
% 1.2145526e-07
% 1.2086969e-07
% 1.2018258e-07
% 1.2080885e-07
% 1.2081281e-07
% 1.2075546e-07
% 1.2119908e-07
% 1.2261300e-07
% 1.2534436e-07
% 1.2640959e-07
% 1.2215129e-07
% 1.0907983e-07
% 8.7132769e-08
% 6.5496849e-08
% 4.8653927e-08
% 3.8614491e-08
% 3.5801506e-08
% 3.2451439e-08
% 3.0297539e-08
% 2.8680056e-08
% 2.7158357e-08
% 2.8170335e-08
% 2.8481233e-08
% 2.9564601e-08
% 3.0449463e-08
% 3.1974338e-08
% 3.6019973e-08
% 4.0316724e-08
% 4.6278682e-08
% 4.9104489e-08
% 4.9653529e-08
% 4.8790913e-08
% 4.7481816e-08
% 4.6222343e-08
% 4.3228642e-08
% 4.3845962e-08
% 4.0778104e-08
% 3.9834198e-08
% 4.1223409e-08
% 3.7333930e-08
% 3.6928320e-08
% 3.6614889e-08
% 3.4177713e-08];
% Bm =[
%     2.86644936177689e-07
% 2.85638999516620e-07
% 2.83628918037259e-07
% 2.81629647996234e-07
% 2.79641149234903e-07
% 2.77663381676735e-07
% 2.75696305327277e-07
% 2.73739880274131e-07
% 2.71794066686914e-07
% 2.69858824817236e-07
% 2.67934114998664e-07
% 2.66019897646696e-07
% 2.64116133258724e-07
% 2.62222782414012e-07
% 2.60339805773657e-07
% 2.58467164080563e-07
% 2.56604818159411e-07
% 2.54752728916624e-07
% 2.52910857340342e-07
% 2.51079164500385e-07
% 2.49257611548225e-07
% 2.47446159716959e-07
% 2.45644770321269e-07
% 2.43853404757399e-07
% 2.42072024503121e-07
% 2.40300591117704e-07
% 2.38539066241882e-07
% 2.36787411597824e-07
% 2.35045588989103e-07
% 2.33313560300664e-07
% 2.31591287498793e-07
% 2.29878732631085e-07
% 2.28175857826414e-07
% 2.26482625294900e-07
% 2.24798997327878e-07
% 2.23124936297868e-07
% 2.21460404658542e-07
% 2.19805364944692e-07
% 2.18159779772198e-07
% 2.16523611837998e-07
% 2.14896823920058e-07
% 2.13279378877333e-07
% 2.11671239649744e-07
% 2.10072369258140e-07
% 2.08482730804268e-07
% 2.06902287470741e-07
% 2.05331002521009e-07
% 2.03768839299320e-07
% 2.02215761230694e-07
% 2.00671731820890e-07
% 1.99136714656370e-07
% 1.97610673404271e-07
% 1.96093571812371e-07
% 1.94585373709056e-07
% 1.93086043003291e-07
% 1.91595543684583e-07
% 1.90113839822950e-07
% 1.88640895568890e-07
% 1.87176675153347e-07
% 1.85721142887681e-07
% 1.84274263163630e-07
% 1.82836000453284e-07
% 1.81406319309045e-07
% 1.79985184363601e-07
% 1.78572560329890e-07
% 1.77168412001066e-07
% 1.75772704250468e-07
% 1.74385402031587e-07
% 1.73006470378031e-07
% 1.71635874403494e-07
% 1.70273579301722e-07
% 1.68919550346482e-07
% 1.67573752891523e-07
% 1.66236152370550e-07
% 1.64906714297185e-07
% 1.63585404264935e-07
% 1.62272187947164e-07
% 1.60967031097049e-07
% 1.59669899547556e-07
% 1.59022006807108e-07
% 1.58376780628436e-07];

load('C:\Users\Owen\OneDrive - Montana State University\Research\Reports\Figures for Disseration\Modeling\BSR.mat','BaSmooth','BmSmooth','HSRLrange');
Ba=BaSmooth;
Bm =BmSmooth;
hsrlRange =HSRLrange;
%%

Ba = fillmissing(interp1(hsrlRange,Ba',rm),'nearest');
Bm = fillmissing(interp1(hsrlRange,Bm',rm),'nearest');

Ba = smoothdata(Ba,'movmean',60);
Bm = smoothdata(Bm,'movmean',60);

% Ba = smoothdata(Ba,'movmean',60);
% Bm = smoothdata(Bm,'movmean',60);




Ba = Ba';
Bm = Bm';

%Ba=zeros(size(Bm));
 %Ba = Bm*1;

 %Bm = zeros(size(Bm));
% 
%  Ba(1:268)=Bm(1:268)*2;
%  Ba(269:(269+149)) = Bm(269:(269+149)) .* linspace(2,1,150)';

%--transmission
Sa = 60;
Ta = exp(-dr*cumtrapz(Ba*Sa));
Sm = (8/3)*pi;
Tm = exp(-dr*cumtrapz(Bm*Sm));

%%
%--spectrum
Spectrum.lambda_online = 769.7958;
Spectrum.lambda_offline = 770.1085;

Spectrum.nu_online = 10^7./Spectrum.lambda_online;                    %[cm-1] Online wavenumber
Spectrum.nu_offline = 10^7./Spectrum.lambda_offline;                  %[cm-1] Offline wavenumber

nuMin = Spectrum.nu_online-0.334*5;                                 %[cm-1] Scan lower bound
nuMax = Spectrum.nu_online+0.334*5;                                 %[cm-1] Scan upper bound
nuMin = Spectrum.nu_online-0.334;                                 %[cm-1] Scan lower bound
nuMax = Spectrum.nu_online+0.334;  
Spectrum.nuBin = 0.00222;                                    %[cm-1] Scan increment
%Spectrum.nuBin = 0.00222/3;  
%Spectrum.nuBin = 4.6667e-4;
nu_scan = (nuMin:Spectrum.nuBin:nuMax);                      %[cm-1](1 x nu) Scan vector

nuMin_off = Spectrum.nu_offline-0.334*5;                                 %[cm-1] Scan lower bound
nuMax_off = Spectrum.nu_offline+0.334*5;                                 %[cm-1] Scan upper bound
nuMin_off = Spectrum.nu_offline-0.334;                                 %[cm-1] Scan lower bound
nuMax_off = Spectrum.nu_offline+0.334; 
nu_scan_off = (nuMin_off:Spectrum.nuBin:nuMax_off);

Spectrum.nu_scan_3D_short = permute(nu_scan, [3 1 2]);       %[cm-1] putting scan in third dimension
Spectrum.nu_scan_3D_short_off = permute(nu_scan_off, [3 1 2]);       %[cm-1] putting scan in third dimension
Spectrum.lambda_scan_3D_short_off = 10^7./Spectrum.nu_scan_3D_short_off;
Spectrum.lambda_scan_3D_short = 10^7./Spectrum.nu_scan_3D_short;

[~,online_index] = min(abs(Spectrum.nu_online-Spectrum.nu_scan_3D_short));
[~,offline_index] = min(abs(Spectrum.nu_offline-Spectrum.nu_scan_3D_short_off));

Spectrum.nu_online = Spectrum.nu_scan_3D_short(:,:,online_index);
Spectrum.nu_offline = Spectrum.nu_scan_3D_short_off(:,:,offline_index);

lambda_online = 10^7./Spectrum.nu_online;
lambda_offline = 10^7./Spectrum.nu_offline;
Spectrum.lambda_online = lambda_online;
Spectrum.lambda_offline = lambda_offline;

%--O2absorption
%%%o2absorption = absorption_O2_770_model_wavenumber(T,P,Spectrum.nu_scan_3D_short,WV);
%o2absorption = absorption_O2_770_PCA(T,P,Spectrum.nu_scan_3D_short,WV);

for i=1:size(Spectrum.nu_scan_3D_short,3)
    o2absorption(:,:,i) = absorption_O2_770_model(T,P,Spectrum.nu_scan_3D_short(:,:,i),WV);
end

%o2absorption_off = absorption_O2_770_model_wavenumber(T,P,Spectrum.nu_scan_3D_short_off ,WV);
%o2absorption_off = absorption_O2_770_PCA_off(T,P,Spectrum.nu_scan_3D_short_off,WV);
o2absorption_off = absorption_O2_770_model(T,P,Spectrum.nu_offline,WV);
absorption = absorption_O2_770_model(T,P,Spectrum.nu_online,WV);
TO2 = exp(-dr*cumtrapz(absorption));
TO2_off = exp(-dr*cumtrapz(o2absorption_off));

absorptionT = -diff(log(TO2))./dr;

absorptionT = zeros(size(TO2));
absorptionT(1,:,:) = (-log(TO2(2)) +log(TO2(1)))./dr;
for iii = 2:size(TO2,1)-1
    absorptionT(iii,:,:) = (-log(TO2(iii+1)) +log(TO2(iii-1)))./dr./2;
end
absorptionT(end,:,:) = (-log(TO2(end)) +log(TO2(end-1)))./dr;

TO2nu = exp(-dr*cumtrapz(o2absorption));
TO2nu_off = exp(-dr*cumtrapz(o2absorption_off));

%%

%--molecularbackscatter

    cB = 1.2;%Brullouion correction to doppler gaussian half width
    cB = -0.01.*((1500+rm')/1000) + 1.2;
    m_air = 28.97/1000./6.02214e23;
    kb = 1.3806e-23;
    c=3e8;

    c_doppler_O2 = m_air*c^2./(8*(Spectrum.nu_online(1)*100).^2*kb);                   %[m^2 K] Doppler coefficient
    doppler_O2_un_ret = ((c_doppler_O2./T/pi).^0.5).*exp(-c_doppler_O2.*(Spectrum.nu_online(1)*100-Spectrum.nu_scan_3D_short*100).^2./T./cB.^2); %[m] Doppler broadended lineshape         

    norm_O2_ret = trapz(doppler_O2_un_ret,3).*Spectrum.nuBin*100;                   %[none] Lineshape integral
    doppler_O2_ret = doppler_O2_un_ret./norm_O2_ret;                       %[m] Normalized doppler lineshape

    c_doppler_O2 = m_air*c^2./(8*(Spectrum.nu_offline(1)*100).^2*kb);                   %[m^2 K] Doppler coefficient
    doppler_O2_un_ret = ((c_doppler_O2./T/pi).^0.5).*exp(-c_doppler_O2.*(Spectrum.nu_offline(1)*100-Spectrum.nu_scan_3D_short_off*100).^2./T./cB.^2); %[m] Doppler broadended lineshape         

    norm_O2_ret = trapz(doppler_O2_un_ret,3).*Spectrum.nuBin*100;                   %[none] Lineshape integral
    doppler_O2_ret_off = doppler_O2_un_ret./norm_O2_ret;                       %[m] Normalized doppler lineshape

gaerosol = zeros(size(Spectrum.nu_scan_3D_short));
gaerosol(online_index) = 1./(Spectrum.nuBin*100);

gaerosol_off = zeros(size(Spectrum.nu_scan_3D_short_off));
gaerosol_off(offline_index) = 1./(Spectrum.nuBin*100);

[Spectrum] = PCAconstrunctionRB2(Spectrum);
[sponts6] = RB_O2_770_PCA(T,P,Spectrum.nu_scan_3D_short,Spectrum);
[sponts6_off] = RB_O2_770_PCA(T,P,Spectrum.nu_scan_3D_short_off,Spectrum);
%[sponts6_off] = RB_O2_770_PCA_offline(T,P,Spectrum.nu_scan_3D_short_off,Spectrum);
doppler_O2_ret = sponts6;
doppler_O2_ret_off = sponts6_off;


% gaerosol = zeros(size(Spectrum.nu_scan_3D_short));
% gaerosol(online_index:online_index+2) = 1./(Spectrum.nuBin*100)/3;
% 
% gaerosol_off = zeros(size(Spectrum.nu_scan_3D_short_off));
% gaerosol_off(offline_index:offline_index+2) = 1./(Spectrum.nuBin*100)/3;

transmitPhotons_on = E_pulse_on./(h*c./(lambda_online*10^-9)) * (pulse_rate/2) * avg_time * 2;
transmitPhotons_off = E_pulse_off./(h*c./(Spectrum.lambda_offline*10^-9)) * (pulse_rate/2) * avg_time* 2;
%%
%%
%spectral purity
% % spectralpurity=1;
% % aseWidth = 10e9;%hz
% % 
% % aseWidth = 100e6;%hz
% % %aseWidth = 0;
% % ase=normpdf(Spectrum.nu_scan_3D_short,Spectrum.nu_online,aseWidth./c./100);
% % ase = ase./trapz(ase)./(Spectrum.nuBin*100);
% % ase_off=normpdf(Spectrum.nu_scan_3D_short_off,Spectrum.nu_offline,aseWidth./c./100);
% % ase_off = ase_off./trapz(ase_off)./(Spectrum.nuBin*100);
% % 
% % laserWidth = gaerosol.*spectralpurity + ase.*(1-spectralpurity);
% % laserWidth_off = gaerosol_off.*spectralpurity + ase_off.*(1-spectralpurity);

spectralpurity =.999;
spectralpurity =.98;
spectralpurity =1;
gaerosol(online_index) = gaerosol(online_index).*spectralpurity;
gaerosol(1:online_index-1) = gaerosol(online_index).*(1-spectralpurity)./(length(gaerosol)-1);
gaerosol(online_index+1:end) = gaerosol(online_index).*(1-spectralpurity)./(length(gaerosol)-1);

gaerosol_off(offline_index) = gaerosol_off(offline_index).*spectralpurity;
gaerosol_off(1:offline_index-1) = gaerosol_off(offline_index).*(1-spectralpurity)./(length(gaerosol)-1);
gaerosol_off(offline_index+1:end) = gaerosol_off(offline_index).*(1-spectralpurity)./(length(gaerosol)-1);

laserWidth = gaerosol;

laserWidth_off = gaerosol_off;

aseWidth = ones(size(laserWidth))./(Spectrum.nuBin*100)./size(laserWidth,3 );

%%
%Overlap

nsPerBin = double(250);                             %[ns] Nanoseconds per bin
NBins = double(560);                                %[none] Number of range bins
rangeBin = (c * nsPerBin(1)*10^-9)/2;               %[m] Create range bin size from speed of light and bin time over 2
rangeMin = -150; %[m]
rangeMin = -300; %[m]
rangeMin = -rangeBin; %[m]
rangeMin = -(c * (1*10^-6))/2;
rangeMin = 0;
rm_raw_o2 = rangeMin:rangeBin:NBins(1)*rangeBin+rangeMin-rangeBin;    %[m] Create range vector
rm_raw_o2 = rm_raw_o2(:);                           %[m] Convert range vector to column vector
r_max = 6000;   %[m] Max range 
%r_max = 10000;   %[m] Max range 
%%
rm_over = rm_raw_o2(rm_raw_o2<=r_max & rm_raw_o2>0);     %[m] Shorten range vector to max range
%load('D:\OneDrive - Montana State University\Research\ARM data analysis\overlap7_21_20.mat','overlap')
% load('C:\Users\Owen\OneDrive - Montana State University\Research\ARM data analysis\overlap7_21_20.mat','overlap')

%eta_O = interp1(rm_over,overlap,rm);
%eta_O = fillmissing(eta_O,'nearest');

  load('OverlapSim_5_31_23.mat','R','OVF','OVF_near')
  eta_O = interp1(R,OVF,rm);
  eta_O = fillmissing(eta_O,'nearest');

eta_O_near = interp1(R,OVF_near,rm);
 eta_O_near = fillmissing(eta_O_near,'nearest');


   load('OverlapSimASE_7_12_23.mat','R','OVF','OVF_near')
   %load('OverlapSimASE_7_35_23.mat','R','OVF','OVF_near')
  eta_OASE = interp1(R,OVF,rm);
  eta_OASE = fillmissing(eta_OASE,'nearest');

eta_O_nearASE = interp1(R,OVF_near,rm);
 eta_O_nearASE = fillmissing(eta_O_nearASE,'nearest');


 %%
 load('C:\Users\Owen\OneDrive - Montana State University\Research\O2 DIAL\analysis\Overlap3_1_24.mat','OVF','OVF_near','R')

 eta_O = fillmissing(interp1(R,OVF(2,:),rm),'nearest');
 eta_O_near = fillmissing(interp1(R,OVF_near(2,:),rm),'nearest');
 spectralpurityO = .98;
  spectralpurityO = .5;

gaerosolASE = ones(size(gaerosol_off)).*gaerosol(online_index).*(1-spectralpurityO)./(length(gaerosol)-1);
gaerosolASE_off = ones(size(gaerosol_off)).*gaerosol_off(offline_index).*(1-spectralpurityO)./(length(gaerosol)-1);

 load('C:\Users\Owen\OneDrive - Montana State University\Research\O2 DIAL\analysis\Overlap23_1_24.mat','OVF2','OVF2_near','R')
 load('C:\Users\Owen\OneDrive - Montana State University\Research\O2 DIAL\analysis\Overlap103_1_24.mat','OVF10','OVF10_near','R')
  % OVF2 =OVF10;
  % OVF2_near = OVF10_near;
 eta_O2 = fillmissing(interp1(R,OVF2(2,:),rm),'nearest');
 eta_O2_near = fillmissing(interp1(R,OVF2_near(2,:),rm),'nearest');


%  load('C:\Users\Owen\OneDrive - Montana State University\Research\Reports\Figures for Disseration\Modeling\Overlap\lens 0607 difflection 9_18.mat','OVF','OVF_near','R')
%  eta_O = fillmissing(interp1(R,OVF(2,:),rm),'nearest');
%  eta_O_near = fillmissing(interp1(R,OVF_near(2,:),rm),'nearest');
% 
%   eta_O(eta_O<=0) = .00001;
% eta_O_near(eta_O_near<=0) =.0001;
% 
%   load('C:\Users\Owen\OneDrive - Montana State University\Research\Reports\Figures for Disseration\Modeling\Overlap\lens 0607 difflection 2 theta 9_18.mat','OVF','OVF_near','R')
%  eta_O2 = fillmissing(interp1(R,OVF(2,:),rm),'nearest');
%  eta_O2_near = fillmissing(interp1(R,OVF_near(2,:),rm),'nearest');
%  eta_O2(eta_O2<=0) = .00001;
% eta_O2_near(eta_O2_near<=0) =.00001;

 %%

Fon  = eta_O'*.98.*Tm.^2.*Ta.^2.*TO2.^2.*(Bm+Ba)./rm'.^2;
Foff  = eta_O'*.98.*Tm.^2.*Ta.^2.*(Bm+Ba)./rm'.^2;

Fon  = eta_O'*.98.*Tm.^2.*Ta.^2.*TO2.^2.*(Bm+Ba)./rm'.^2 + eta_O_near'*.02.*Tm.^2.*Ta.^2.*TO2.^2.*(Bm+Ba)./rm'.^2;
Foff  = eta_O'*.98.*Tm.^2.*Ta.^2.*(Bm+Ba)./rm'.^2 + eta_O_near'*.02.*Tm.^2.*Ta.^2.*(Bm+Ba)./rm'.^2;

Fon  = eta_O'*.98.*Tm.^2.*Ta.^2.*TO2.^2.*(Bm+Ba)./rm'.^2 + eta_O_near'*.02.*Tm.^2.*Ta.^2.*TO2.^2.*(Bm+Ba)./rm'.^2;
Foff  = eta_O'*.98.*Tm.^2.*Ta.^2.*TO2_off.^2.*(Bm+Ba)./rm'.^2 + eta_O_near'*.02.*Tm.^2.*Ta.^2.*TO2_off.^2.*(Bm+Ba)./rm'.^2;

load('CalibrationData\TransmissionData20220809.mat','Data_Wavelength')
eta_T_on = interp1(Data_Wavelength.lambda_on.*10^9,Data_Wavelength.Tc_on./max(Data_Wavelength.Tc_on),Spectrum.lambda_scan_3D_short);
eta_T_off = interp1(Data_Wavelength.lambda_off.*10^9,Data_Wavelength.Tc_off./max(Data_Wavelength.Tc_off),Spectrum.lambda_scan_3D_short_off);


FonSpectrumUp = trapz(Tm.*Ta.*TO2nu.*laserWidth,3).*Spectrum.nuBin*100;
FonSpectrumDown = trapz(eta_T_on.*Tm.*Ta.*TO2nu.*(Bm.*(convn(doppler_O2_ret,laserWidth,'same').*Spectrum.nuBin.*100) ...
    + Ba.*laserWidth),3).*Spectrum.nuBin*100;

FoffSpectrumUp = trapz(Tm.*Ta.*TO2nu_off.*laserWidth_off,3).*Spectrum.nuBin*100;
FoffSpectrumDown = trapz(eta_T_off.*Tm.*Ta.*TO2nu_off.*(Bm.*(convn(doppler_O2_ret_off,laserWidth_off,'same').*Spectrum.nuBin.*100) ...
    + Ba.*laserWidth_off),3).*Spectrum.nuBin*100;

FonSpectrumUpASE = trapz(Tm.*Ta.*TO2nu.*aseWidth,3).*Spectrum.nuBin*100;
FonSpectrumDownASE = trapz(eta_T_on.*Tm.*Ta.*TO2nu.*(Bm.*(convn(doppler_O2_ret,aseWidth,'same').*Spectrum.nuBin.*100) ...
    + Ba.*aseWidth),3).*Spectrum.nuBin*100;


FoffSpectrumUpASE = trapz(Tm.*Ta.*TO2nu_off.*aseWidth,3).*Spectrum.nuBin*100;
FoffSpectrumDownASE = trapz(eta_T_off.*Tm.*Ta.*TO2nu_off.*(Bm.*(convn(doppler_O2_ret_off,aseWidth,'same').*Spectrum.nuBin.*100) ...
    + Ba.*aseWidth),3).*Spectrum.nuBin*100;


FonSpectrumUpASEO = trapz(Tm.*Ta.*TO2nu.*gaerosolASE,3).*Spectrum.nuBin*100;
FonSpectrumDownASEO = trapz(eta_T_on.*Tm.*Ta.*TO2nu.*(Bm.*(convn(doppler_O2_ret,gaerosolASE,'same').*Spectrum.nuBin.*100) ...
    + Ba.*laserWidth),3).*Spectrum.nuBin*100;

FoffSpectrumUpASEO = trapz(Tm.*Ta.*TO2nu_off.*gaerosolASE_off,3).*Spectrum.nuBin*100;
FoffSpectrumDownASEO = trapz(eta_T_off.*Tm.*Ta.*TO2nu_off.*(Bm.*(convn(doppler_O2_ret_off,gaerosolASE_off,'same').*Spectrum.nuBin.*100) ...
    + Ba.*laserWidth_off),3).*Spectrum.nuBin*100;

%%
%Afterpulse

x = rm/(250e-9*3e8/2);

a = 148.0197;
b = -0.6118;
c = 13.7442;
d = -0.0936;

% a=2.1315e7;
% b=-0.6118;
% c=1.9792e6;
% d=-0.0936;
Correction_NC_off = a*exp(b*(x-5))+c*exp(d*(x-5));

a = 152.7929;
b = -0.6431;
c = 16.1950;
d = -0.1081;

% a=2.2002e7;
% b=-0.6431;
% c=2.3321e6;
% d=-0.1081;

Correction_NC_on= a*exp(b*(x-5))+c*exp(d*(x-5));

% Correction_NC_off = Correction_NC_off'*avg_time * 2*2;
% Correction_NC_on  = Correction_NC_on'*avg_time * 2*2;

Correction_NC_off = Correction_NC_off'*avg_time * 2*2/1000;
Correction_NC_on  = Correction_NC_on'*avg_time * 2*2/1000;

%%
N_on  = transmitPhotons_on  .* eta_R .* eta_D .* A .* dr .* Fon *2.5*2;
N_off = transmitPhotons_on .* eta_R .* eta_D .* A .* dr .* Foff *2.5*2;

load('Afterpulsing_correction_09092022.mat','Correction_Nc_on','Correction_Nc_off')
correctionRange = 0:(250e-9*3e8/2):(length(Correction_Nc_on)-1)*(250e-9*3e8/2);
Correction_Nc_on = interp1(correctionRange,Correction_Nc_on,rm)';
Correction_Nc_off = interp1(correctionRange,Correction_Nc_off,rm)';

%Correction_Nc_on = [Correction_Nc_on'; mean(Correction_Nc_on(:,end-20:end)).*ones(length(N_on)-length(Correction_Nc_on),1)];
%Correction_Nc_off = [Correction_Nc_off'; mean(Correction_Nc_off(:,end-20:end)).*ones(length(N_off)-length(Correction_Nc_off),1)];
N_onAfterpulse = N_on + movmean((Correction_Nc_on),5)-mean(Correction_Nc_on(end-20*4,end));
N_offAfterpulse = N_off + movmean((Correction_Nc_off),5)-mean(Correction_Nc_off(end-20*4,end));

N_onSpectrum = eta_O'.*.98.*FonSpectrumUp.*FonSpectrumDown.*eta_R.*eta_D.*A.*dr.*transmitPhotons_on./rm'.^2 + eta_O_near'.*.02.*FonSpectrumUp.*FonSpectrumDown.*eta_R.*eta_D.*A.*dr.*transmitPhotons_on./rm'.^2;
N_offSpectrum = eta_O'.*.98.*FoffSpectrumUp.*FoffSpectrumDown.*eta_R.*eta_D.*A.*dr.*transmitPhotons_off./rm'.^2 + eta_O_near'.*.02.*FoffSpectrumUp.*FoffSpectrumDown.*eta_R.*eta_D.*A.*dr.*transmitPhotons_off./rm'.^2;


N_onSpectrumASE = eta_OASE'.*.98.*FonSpectrumUpASE.*FonSpectrumDownASE.*eta_R.*eta_D.*A.*dr.*transmitPhotons_on*.01./rm'.^2 + eta_O_nearASE'.*.02.*FonSpectrumUpASE.*FonSpectrumDownASE.*eta_R.*eta_D.*A.*dr.*transmitPhotons_on*.01./rm'.^2;
N_offSpectrumASE = eta_OASE'.*.98.*FoffSpectrumUpASE.*FoffSpectrumDownASE.*eta_R.*eta_D.*A.*dr.*transmitPhotons_off*.01./rm'.^2 + eta_O_nearASE'.*.02.*FoffSpectrumUpASE.*FoffSpectrumDownASE.*eta_R.*eta_D.*A.*dr.*transmitPhotons_off*.01./rm'.^2;


N_onSpectrumASEO = eta_O2'.*.98.*FonSpectrumUpASEO.*FonSpectrumDownASEO.*eta_R.*eta_D.*A.*dr.*transmitPhotons_on*.01./rm'.^2 + eta_O2_near'.*.02.*FonSpectrumUpASEO.*FonSpectrumDownASEO.*eta_R.*eta_D.*A.*dr.*transmitPhotons_on*.01./rm'.^2;
N_offSpectrumASEO = eta_O2'.*.98.*FoffSpectrumUpASEO.*FoffSpectrumDownASEO.*eta_R.*eta_D.*A.*dr.*transmitPhotons_off*.01./rm'.^2 + eta_O2_near'.*.02.*FoffSpectrumUpASEO.*FoffSpectrumDownASEO.*eta_R.*eta_D.*A.*dr.*transmitPhotons_off*.01./rm'.^2;


N_onSpectrumAfterpulse = N_onSpectrum+Correction_NC_on;
N_offSpectrumAfterpulse = N_offSpectrum+Correction_NC_off;
%%

BSR = (Ba+Bm)./Bm;
smoothingLength = 300;%pulse length
smoothingLength = 150;%pulse length
smoothPoints = round(smoothingLength/dr);
%smoothPoints=1;
smoothVector = ones(smoothPoints,1)/smoothPoints;


%%% Summing over pulse length
for iii = 1:length(rm)
    if iii <=smoothPoints
        
        N_onPulse(iii,:) = mean(N_on(1:iii,:),1);
        N_offPulse(iii,:) = mean(N_off(1:iii,:),1);
    
        N_onSpectrumPulse(iii,:) = mean(N_onSpectrum(1:iii,:),1);
        N_offSpectrumPulse(iii,:) = mean(N_offSpectrum(1:iii,:),1);

        BSRshift(iii,:) = mean(BSR(1:iii,:),1);


    else
        N_onPulse(iii,:) = sum(N_on(iii-smoothPoints:iii,:),1)./smoothPoints;
        N_offPulse(iii,:) = sum(N_off(iii-smoothPoints:iii,:),1)./smoothPoints;
    
        N_onSpectrumPulse(iii,:) = sum(N_onSpectrum(iii-smoothPoints:iii,:),1)./smoothPoints;
        N_offSpectrumPulse(iii,:) = sum(N_offSpectrum(iii-smoothPoints:iii,:),1)./smoothPoints;

        BSRshift(iii,:) = sum(BSR(iii-smoothPoints:iii,:),1)./smoothPoints;
        
    end
end

N_onPulseAfterpulse = N_onPulse- movmean((Correction_Nc_on),5)-mean(Correction_Nc_on(end-20*4,end));
N_offPulseAfterpulse = N_offPulse- movmean((Correction_Nc_off),5)-mean(Correction_Nc_off(end-20*4,end));


% BSRshift = [ones(smoothPoints,1).*BSR(1); BSR(1:(end-smoothPoints))];

N_onPulse = conv(N_on,smoothVector,'same');
N_offPulse = conv(N_off,smoothVector,'same');

% N_onPulse(1:end-smoothPoints./2) = N_onPulse(smoothPoints./2+1:end);
% N_offPulse(1:end-smoothPoints./2) = N_offPulse(smoothPoints./2+1:end);

% N_onPulse(smoothPoints./2:end)=N_onPulse(1:end-smoothPoints./2+1) ;
% N_offPulse(smoothPoints./2:end) = N_offPulse(1:end-smoothPoints./2+1) ;

 % N_onPulse(2:end)=N_onPulse(1:end-1) ;
 % N_offPulse(2:end) = N_offPulse(1:end-1) ;

% N_onPulse(1:end-smoothPoints) = N_onPulse(smoothPoints+1:end);
% N_offPulse(1:end-smoothPoints) = N_offPulse(smoothPoints+1:end);

% N_onPulse(1:end-1) = N_onPulse(2:end);
% N_offPulse(1:end-1) = N_offPulse(2:end);

N_onSpectrumPulse = conv(N_onSpectrum,smoothVector,'same');
N_offSpectrumPulse = conv(N_offSpectrum,smoothVector,'same');

% N_onSpectrumPulse(1:end-smoothPoints./2) = N_onSpectrumPulse(smoothPoints./2+1:end);
% N_offSpectrumPulse(1:end-smoothPoints./2) = N_offSpectrumPulse(smoothPoints./2+1:end);

 N_onSpectrumPulse(2:end)=N_onSpectrumPulse(1:end-1) ;
 N_offSpectrumPulse(2:end) = N_offSpectrumPulse(1:end-1) ;

 N_onSpectrumPulseASEO = conv(N_onSpectrum+N_onSpectrumASEO,smoothVector,'same');
N_offSpectrumPulseASEO = conv(N_offSpectrum+N_offSpectrumASEO,smoothVector,'same');
 N_onSpectrumPulseASEO = conv(N_onSpectrum+N_onSpectrumASEO+N_offSpectrumASEO,smoothVector,'same');
N_offSpectrumPulseASEO = conv(N_offSpectrum+N_offSpectrumASEO+N_onSpectrumASEO,smoothVector,'same');

BSRPulse = conv(BSR,smoothVector,'same');

%BSRPulse(1:end-smoothPoints./2) = BSR(smoothPoints./2+1:end);

%%
%Background correction
%OnlineBackground nighttime 2021 
% % N_onSpectrumPulse= N_onSpectrumPulse-2021*.01;
% % N_offSpectrumPulse= N_offSpectrumPulse+2021*.01;
% % %%
% % N_onSpectrum= N_onSpectrum-2021*.01;
% % N_offSpectrum= N_offSpectrum+2021*.01;

%%
[a_0] = alpha_0(N_on,N_off,dr);
[a_0Pulse] = alpha_0(N_onPulse,N_offPulse,dr);
[a_0Spectrum] = alpha_0(N_onSpectrum,N_offSpectrum,dr);
[a_0SpectrumPulse] = alpha_0(N_onSpectrumPulse,N_offSpectrumPulse,dr);
[a_0Afterpulse] = alpha_0(N_onAfterpulse,N_offAfterpulse,dr);
[a_0PulseAfterpulse] = alpha_0(N_onPulseAfterpulse,N_offPulseAfterpulse,dr);

[a_0SpectrumAfterpulse] = alpha_0(N_onSpectrumAfterpulse,N_offSpectrumAfterpulse,dr);
[a_0SpectrumPulseAfterpulse] = alpha_0(N_onSpectrumPulse+Correction_NC_on,N_offSpectrumPulse+Correction_NC_off,dr);

[a_0SpectrumASE] = alpha_0(N_onSpectrum+N_onSpectrumASEO,N_offSpectrum+N_offSpectrumASEO,dr);
[a_0SpectrumASE] = alpha_0(N_onSpectrum+N_onSpectrumASEO+N_offSpectrumASEO,N_offSpectrum+N_offSpectrumASEO+N_onSpectrumASEO,dr);

[a_0SpectrumPulseASE] = alpha_0(N_onSpectrumPulseASEO,N_offSpectrumPulseASEO,dr);


a_0 = a_0+o2absorption_off;
a_0Pulse = a_0Pulse+o2absorption_off;
a_0Spectrum = a_0Spectrum+o2absorption_off;
a_0SpectrumPulse = a_0SpectrumPulse+o2absorption_off;
a_0Afterpulse = a_0Afterpulse+o2absorption_off;
a_0PulseAfterpulse = a_0PulseAfterpulse+o2absorption_off;
a_0Spectrumfterpulse = a_0SpectrumAfterpulse+o2absorption_off;
% 
% figure(1)
% plot(a_0,rm)
% hold on
% plot(a_0Pulse,rm)
% plot(a_0Spectrum,rm)
% plot(a_0SpectrumPulse,rm)
% plot(a_0Afterpulse,rm)
% plot(a_0PulseAfterpulse,rm)
% plot(o2absorption(:,:,online_index),rm,'--')
% hold off
% legend('a_0','a_0Pulse','a_0Spectrum','a_0SpectrumPulse','a_0 afterpulse','absorption')

figure(2)
plot(a_0-o2absorption(:,:,online_index),rm)
hold on
plot(a_0Pulse-o2absorption(:,:,online_index),rm,'.-')
plot(a_0Spectrum-o2absorption(:,:,online_index),rm,'--')
plot(a_0SpectrumPulse-o2absorption(:,:,online_index),rm,'--')

plot(a_0SpectrumASE-o2absorption(:,:,online_index),rm,'.')
%plot(a_0Afterpulse-o2absorption(:,:,online_index),rm)
%plot(a_0PulseAfterpulse-o2absorption(:,:,online_index),rm,'--')
hold off
legend('a_0','a_0Pulse','a_0Spectrum','a_0SpectrumPulse')

figure(22)
plot((a_0-o2absorption(:,:,online_index))./o2absorption(:,:,online_index)*100,rm)
hold on
plot((a_0Pulse-o2absorption(:,:,online_index))./o2absorption(:,:,online_index)*100,rm,'.-')
plot((a_0Spectrum-o2absorption(:,:,online_index))./o2absorption(:,:,online_index)*100,rm,'--')
plot((a_0SpectrumPulse-o2absorption(:,:,online_index))./o2absorption(:,:,online_index)*100,rm,'--')
%plot(a_0Afterpulse-o2absorption(:,:,online_index),rm)
%plot(a_0PulseAfterpulse-o2absorption(:,:,online_index),rm,'--')
hold off
legend('a_0','a_0Pulse','a_0Spectrum','a_0SpectrumPulse')

%%

%--correction
altitude = 1.5719;
Model.T = T;
Model.P = P;
Model.WV = 0;
BSR = (Ba+Bm)./Bm;

Options.oversample=1;
Range.rm = rm';
Range.i_range = length(rm);
Range.rangeBin = rm(2)-rm(1);
Time.ts = 1;
Time.i_time = 1;
Spectrum.i_scan_3D_short = length(Spectrum.nu_scan_3D_short);
Spectrum.online_index = online_index;
Spectrum.offline_index = offline_index;

[alpha_final,alpha_1_rawSpectrum,alpha_2_rawSpectrum,Spectrum] = pertAbsorption(a_0Spectrum, eta_T_on,eta_T_off, Model, Range, Time, Spectrum, BSR, 0,0, Options, 0,false);
Alpha_totalSpectrum = a_0Spectrum+alpha_1_rawSpectrum+alpha_2_rawSpectrum;

[alpha_final,alpha_1_rawSpectrumAfterpulse,alpha_2_rawSpectrumAfterpulse,Spectrum] = pertAbsorption(a_0SpectrumAfterpulse, eta_T_on,eta_T_off, Model, Range, Time, Spectrum, BSR, 0,0, Options, 0,false);
Alpha_totalSpectrumAfterpulse = a_0SpectrumAfterpulse+alpha_1_rawSpectrumAfterpulse+alpha_2_rawSpectrumAfterpulse;

[alpha_final,alpha_1_rawSpectrumPulseAfterpulse,alpha_2_rawSpectrumPulseAfterpulse,Spectrum] = pertAbsorption(a_0SpectrumPulseAfterpulse, eta_T_on,eta_T_off, Model, Range, Time, Spectrum, BSRPulse, 0,0, Options, 0,false);
Alpha_totalSpectrumPulseAfterpulse = a_0SpectrumPulseAfterpulse+alpha_1_rawSpectrumPulseAfterpulse+alpha_2_rawSpectrumPulseAfterpulse;


[alpha_final,alpha_1_rawSpectrumASE,alpha_2_rawSpectrumASE,Spectrum] = pertAbsorption(a_0SpectrumASE, eta_T_on,eta_T_off, Model, Range, Time, Spectrum, BSR, 0,0, Options, 0,false);
Alpha_totalSpectrumASE = a_0SpectrumASE+alpha_1_rawSpectrumASE+alpha_2_rawSpectrumASE;

[alpha_final,alpha_1_rawSpectrumPulseASE,alpha_2_rawSpectrumPulseASE,Spectrum] = pertAbsorption(a_0SpectrumPulseASE, eta_T_on,eta_T_off, Model, Range, Time, Spectrum, BSRPulse, 0,0, Options, 0,false);
Alpha_totalSpectrumPulseASE = a_0SpectrumPulseASE+alpha_1_rawSpectrumPulseASE+alpha_2_rawSpectrumPulseASE;


[alpha_final,alpha_1_rawSpectrumPulse,alpha_2_rawSpectrumPulse,Spectrum] = pertAbsorption(a_0SpectrumPulse, eta_T_on,eta_T_off, Model, Range, Time, Spectrum, BSR, 0,0, Options, 0,false);
Alpha_totalSpectrumPulse = a_0SpectrumPulse+alpha_1_rawSpectrumPulse+alpha_2_rawSpectrumPulse;

[alpha_final,alpha_1_rawSpectrumPulseBSR,alpha_2_rawSpectrumPulseBSR,Spectrum] = pertAbsorption(a_0SpectrumPulse, eta_T_on,eta_T_off, Model, Range, Time, Spectrum, BSRPulse, 0,0, Options, 0,false);
Alpha_totalSpectrumPulseBSRShift = a_0SpectrumPulse+alpha_1_rawSpectrumPulseBSR+alpha_2_rawSpectrumPulseBSR;

figure(3)
plot(a_0,rm)
hold on
plot(a_0Pulse,rm)
plot(Alpha_totalSpectrum,rm)
plot(Alpha_totalSpectrumPulse,rm,'--')
plot(Alpha_totalSpectrumPulseBSRShift,rm,'.-')
plot(o2absorption(:,:,online_index),rm,'--')
hold off
legend('a_t','a_tPulse','a_tSpectrum','a_tSpectrumPulse','bsrShiftSpectrumpulse','absorption')

figure(4)
% plot(a_0-absorption,rm)
% hold on
% plot(a_0Pulse-absorption,rm)
plot((Alpha_totalSpectrum-absorption)./absorption*100,rm/1000,'linewidth',2)
hold on
%plot(Alpha_totalSpectrumPulse-absorption,rm/1000,'.-','linewidth',2)
plot((Alpha_totalSpectrumPulseBSRShift-absorption)./absorption*100,rm/1000,'--','linewidth',2)
%plot(o2absorption(:,:,online_index)-absorption,rm/1000,'--')
plot((Alpha_totalSpectrumASE-absorption)./absorption*100,rm/1000,'-')
plot((Alpha_totalSpectrumPulseASE-absorption)./absorption*100,rm/1000,'--')
%plot((Alpha_totalSpectrumAfterpulse-absorption)./absorption*100,rm/1000,'LineWidth',2)
%plot((Alpha_totalSpectrumPulseAfterpulse-absorption)./absorption*100,rm/1000,'--','LineWidth',2)
hold off
legend('a_t','a_tPulse','a_tSpectrum','a_tSpectrumPulse','bsrShiftSpectrumpulse','absorption')
legend('Absorption','Absorption Pulse','Afterpulse','Afterpulse Pulse')
legend('Absorption','Absorption Pulse','ASE')
grid on
ylim([0 6])
xlim([-5 5])
xlabel('Absorption % Difference (\alpha_{retrieved}-\alpha)')
ylabel('Range (km)')

figure(44)
plot(alpha_1_rawSpectrum,rm)
hold on
plot(alpha_1_rawSpectrumPulse,rm)
%plot(alpha_1_rawSpectrum,rm)
hold off
legend('a_tSpectrum','a_tSpectrumPulse')

%%

[T_final,Lapse,Ts_fit,P_final,mean_lapse_rate,exclusion,Titer] =                        temperatureRetrieval(T,1,rm',0,WV,Spectrum.nu_online,a_0,0,logical(zeros(size(rm'))),T(1),P(1),-6.5/1000);
[T_finalPulse,Lapse,Ts_fit,P_finalPulse,mean_lapse_rate,exclusion,Titer] =                   temperatureRetrieval(T,1,rm',0,WV,Spectrum.nu_online,a_0Pulse,0,logical(zeros(size(rm'))),T(1),P(1),-6.5/1000);
[T_finalSpectrum,Lapse,Ts_fit,P_finalSpectrum,mean_lapse_rate,exclusion,Titer] =                temperatureRetrieval(T,1,rm',0,WV,Spectrum.nu_online,Alpha_totalSpectrum,0,logical(zeros(size(rm'))),T(1),P(1),-6.5/1000);
[T_finalSpectrumAfterpulse,Lapse,Ts_fit,P_finalSpectrum,mean_lapse_rate,exclusion,Titer] =                temperatureRetrieval(T,1,rm',0,WV,Spectrum.nu_online,Alpha_totalSpectrumAfterpulse,0,logical(zeros(size(rm'))),T(1),P(1),-6.5/1000);
[T_finalSpectrumPulseAfterpulse,Lapse,Ts_fit,P_finalSpectrum,mean_lapse_rate,exclusion,Titer] =                temperatureRetrieval(T,1,rm',0,WV,Spectrum.nu_online,Alpha_totalSpectrumPulseAfterpulse,0,logical(zeros(size(rm'))),T(1),P(1),-6.5/1000);


[T_finalSpectrumPulse,Lapse,Ts_fit,P_finalSpectrumPulse,mean_lapse_rate,exclusion,Titer] =           temperatureRetrieval(T,1,rm',0,WV,Spectrum.nu_online,Alpha_totalSpectrumPulse,0,logical(zeros(size(rm'))),T(1),P(1),-6.5/1000);
[T_finalSpectrumPulseBSRShift,Lapse,Ts_fit,P_finalPulseBSRShift,mean_lapse_rate,exclusion,Titer] =   temperatureRetrieval(T,1,rm',0,WV,Spectrum.nu_online,Alpha_totalSpectrumPulseBSRShift,0,logical(zeros(size(rm'))),T(1),P(1),-6.5/1000);

[T_finalSpectrumASE,Lapse,Ts_fit,P_finalSpectrum,mean_lapse_rate,exclusion,Titer] =                temperatureRetrieval(T,1,rm',0,WV,Spectrum.nu_online,Alpha_totalSpectrumASE,0,logical(zeros(size(rm'))),T(1),P(1),-6.5/1000);
[T_finalSpectrumPulseASE,Lapse,Ts_fit,P_finalSpectrum,mean_lapse_rate,exclusion,Titer] =                temperatureRetrieval(T,1,rm',0,WV,Spectrum.nu_online,Alpha_totalSpectrumPulseASE,0,logical(zeros(size(rm'))),T(1),P(1),-6.5/1000);


figure(5)
plot(T_final,rm)
hold on
plot(T_finalPulse,rm)
plot(T_finalSpectrum,rm)
plot(T_finalSpectrumPulse,rm)
plot(T_finalSpectrumPulseBSRShift,rm)
plot(T,rm,'--')
hold off
legend('Tfinal','TfinalPulse','TfinalSpectrum','TfinalSpectrumPulse','shift','T')

figure(6)
plot(T-T_final,rm)
hold on
plot(T-T_finalPulse,rm)
plot(T-T_finalSpectrum,rm)
plot(T-T_finalSpectrumPulse,rm)
plot(T-T_finalSpectrumPulseBSRShift,rm,'--','linewidth',2)
plot(T-T_finalSpectrumASE,rm,'.-')
plot(T-T_finalSpectrumAfterpulse,rm,'.')
plot(T-T_finalSpectrumPulseAfterpulse,rm,'.')
%plot(T(1:end-smoothPoints)-T_finalSpectrumPulseBSRShift((smoothPoints+1):end),rm(1:end-smoothPoints),'--')
hold off
legend('Tfinal','TfinalPulse','TfinalSpectrum','TfinalSpectrumPulse','BSR shift')

% figure(7)
% plot(T-T_final,rm)
% hold on
% plot(T(1:end-smoothPoints)-T_finalPulse((smoothPoints+1):end),rm(1:end-smoothPoints))
% plot(T-T_finalSpectrum,rm)
% plot(T(1:end-smoothPoints)-T_finalSpectrumPulse((smoothPoints+1):end),rm(1:end-smoothPoints))
% %plot(T(1:end-smoothPoints)-T_finalSpectrumPulseBSRShift((smoothPoints+1):end),rm(1:end-smoothPoints),'--')
% %plot(T(1:end-smoothPoints)-T_finalSpectrumPulseBSRShift((smoothPoints+1):end),rm(1:end-smoothPoints),'--')
% hold off
% legend('Tfinal','TfinalPulse','TfinalSpectrum','TfinalSpectrumPulse')

figure(8)
plot(BSR,rm/1000,'linewidth',2)
xlabel('BSR')
ylabel('Range (km)')
grid on
ylim([0 6])

figure(9)
plot(eta_O,rm/1000,'linewidth',2)
hold on
plot(eta_O_near,rm/1000,'linewidth',2)
% plot(eta_OASE,rm/1000,'--','linewidth',2)
% plot(eta_O_nearASE,rm/1000,'--','linewidth',2)
plot(eta_O2,rm/1000,'--','linewidth',2)
plot(eta_O2_near,rm/1000,'--','linewidth',2)
hold off
xlabel('Lidar overlap function')
ylabel('Range (km)')
legend('Primary telescope','Secondary telescope','Primary telescope ASE','Secondary telescope ASE')
ylim([0 6])

figure(10)
% plot(N_on,rm/1000,'linewidth',2)
% hold on
semilogx(N_onSpectrum,rm/1000,'linewidth',2)
hold on
semilogx(N_offSpectrum,rm/1000,'linewidth',2)
%semilogx(N_onSpectrumPulse,rm/1000,'--','linewidth',2)
%semilogx(N_offSpectrumPulse,rm/1000,'--','linewidth',2)
hold off
xlim([500 0.5e6])
ylim([0 6])
legend('Online','Offline','Online Pulse','Offline Pulse')
xlabel('Modeled Signal Counts')
ylabel('Range (rkm)')
grid on

figure(11)
plot(P_final-P,rm)
hold on
plot(P_finalPulse-P,rm)
plot(P_finalSpectrum-P,rm,'--')
plot(P_finalSpectrumPulse-P,rm,'--')
hold off
legend('P','P_pulse','P spectrum','P Spectrum pulse')


figure(444)
plot(T_finalSpectrum-T,rm/1000,'linewidth',2)
hold on
plot(T_finalSpectrumPulseBSRShift-T,rm/1000,'--','linewidth',2)
% plot(T_finalSpectrumAfterpulse-T,rm/1000,'linewidth',2)
% plot(T_finalSpectrumPulseAfterpulse-T,rm/1000,'--','linewidth',2)
plot(T_finalSpectrumASE-T,rm/1000,'-')
plot(T_finalSpectrumPulseASE-T,rm/1000,'-')
hold off
legend('a_t','a_tPulse','a_tSpectrum','a_tSpectrumPulse','bsrShiftSpectrumpulse','absorption')
legend('Temperature','Temperature Pulse','Afterpulse','Afterpulse Pulse')
grid on
ylim([0 6])
xlim([-3 3])
xlabel('Temperature Difference (T_{retrieved}-T)')
ylabel('Range (km)')



