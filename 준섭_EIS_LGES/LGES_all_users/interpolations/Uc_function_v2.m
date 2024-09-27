function [ U_avg ] = Uc_function_v2( x,chg )
% This function interpolates the equilibrium potential curve of cathode active
% material. Averaging the cuves obtanied in charging and discharging.

%% Discharge data: "NMC_22FKV08_OCV_HC2211_Dchg.txt"
data_dis = [0.323954472	4.347948074
            0.330265803	4.332117081
            0.33657736	4.316770554
            0.342888667	4.30190897
            0.349199981	4.287531853
            0.355511308	4.273155212
            0.361822544	4.259100914
            0.368133836	4.24520874
            0.374445295	4.230993271
            0.380756713	4.217585564
            0.387067964	4.203854561
            0.393379498	4.190123558
            0.399690765	4.176392555
            0.406002021	4.162984848
            0.412313298	4.149253845
            0.418624612	4.136007786
            0.424936155	4.122600079
            0.431247668	4.109354019
            0.437559087	4.096430779
            0.443870416	4.083507538
            0.450181847	4.070422649
            0.456493186	4.057984352
            0.462804428	4.045707226
            0.469115781	4.032945633
            0.475427087	4.020991802
            0.481738531	4.008553028
            0.488050064	3.996760607
            0.494361601	3.984968185
            0.500672894	3.973337412
            0.506984142	3.961544991
            0.51329546	3.950237274
            0.519606813	3.93876791
            0.525918074	3.92794466
            0.532229607	3.916636944
            0.53854113	3.905652285
            0.544852568	3.89515233
            0.551164117	3.884975195
            0.557475353	3.875282764
            0.563786839	3.866075039
            0.57009835	3.857028961
            0.57640977	3.848305702
            0.582721077	3.840067148
            0.589032465	3.832636356
            0.595343956	3.825043917
            0.601655186	3.818420887
            0.607966658	3.811797857
            0.614278078	3.805659294
            0.620589443	3.800166845
            0.626900931	3.794674635
            0.633212187	3.7898283
            0.639523622	3.784982204
            0.645835087	3.780459166
            0.652146368	3.776582003
            0.658457868	3.772705078
            0.66476932	3.768989801
            0.671080566	3.765597343
            0.677391814	3.762366533
            0.683703227	3.759135723
            0.690014531	3.756066561
            0.696325802	3.753320456
            0.702637137	3.750089645
            0.70894866	3.747828007
            0.715259969	3.745081902
            0.721571461	3.742335796
            0.727882992	3.73991251
            0.73419454	3.737328053
            0.740505769	3.735066414
            0.746817148	3.732481718
            0.753128689	3.729897261
            0.75944014	3.727312565
            0.76575139	3.724889517
            0.772062645	3.722304821
            0.778373872	3.719881773
            0.784685225	3.717135429
            0.790996659	3.714550972
            0.797307945	3.711966276
            0.803619182	3.709058523
            0.809930731	3.706312418
            0.816242222	3.703404665
            0.82255372	3.70065856
            0.828865154	3.697589159
            0.835176646	3.694519997
            0.841487879	3.691289186
            0.847799105	3.688058376
            0.854110661	3.684504509
            0.860421912	3.681112289
            0.866733291	3.677558422
            0.873044684	3.673681259
            0.879356238	3.669642925
            0.885667849	3.665604353
            0.891979377	3.661081314
            0.898290627	3.656558037
            0.904602135	3.651550531
            0.910913505	3.646542788
            0.917225001	3.641050339
            0.923536403	3.635073423
            0.92984784	3.628450155
            0.936159193	3.621342421
            0.942470549	3.613104105
            0.948781822	3.604380846
            0.95509307	3.594688416
            0.961404313	3.584672928
            0.96771562	3.570296049
            0.97402686	3.554303646
            0.980338087	3.527649403
            0.985156247	3.49744153
            0.988226584	3.463356733
            0.99022484	3.435894966
            0.991811782	3.400840998
            0.993119486	3.365140676
            0.994241455	3.316355705
            0.9951849	3.247378349
            0.995949133	3.19197011
            0.996589755	3.125738859
            0.997105738	3.109746456
            0.997581779	3.096015692
            0.998026181	3.089069366
            0.998454926	3.07776165
            0.998861784	3.067261457
            0.999250439	3.062092304
            0.99962974	3.059669256
            1	3.046907425];

%% Data charging: "NMC_22FKV08_OCV_HC2211_Chg.txt"

data_cha = [1	3.504226208
            0.993621644	3.562865257
            0.98724315	3.583542347
            0.980864705	3.598403931
            0.974486303	3.608742476
            0.968107958	3.617465496
            0.961729645	3.625219345
            0.955351281	3.632004023
            0.948972948	3.637981176
            0.942594645	3.643958092
            0.936216267	3.649288893
            0.929837668	3.654296637
            0.923459348	3.65930438
            0.917080774	3.66366601
            0.910702145	3.668350458
            0.904323549	3.67238903
            0.897944956	3.676589012
            0.891566344	3.680465937
            0.885187619	3.684342861
            0.878808869	3.687896729
            0.872430137	3.691612244
            0.866051788	3.695004702
            0.859673403	3.698235273
            0.853295058	3.701304674
            0.846916693	3.704696894
            0.840538165	3.707604647
            0.834159625	3.710350752
            0.827781087	3.713258505
            0.821402702	3.716004848
            0.815023988	3.718912363
            0.808645377	3.721335649
            0.802266825	3.723920107
            0.795888183	3.726666451
            0.789509582	3.729089499
            0.783131149	3.732158661
            0.776752635	3.734581709
            0.770374029	3.737004995
            0.763995697	3.739589453
            0.757617281	3.742012739
            0.751238607	3.744435787
            0.744859882	3.747181892
            0.738481568	3.749766588
            0.732103254	3.752189636
            0.72572493	3.754935741
            0.719346329	3.757358789
            0.712967689	3.760266542
            0.706589005	3.763174295
            0.700210635	3.766082048
            0.693832297	3.769312859
            0.687453986	3.772543669
            0.681075237	3.775935888
            0.674696889	3.779651403
            0.668318544	3.783528328
            0.661940171	3.78772831
            0.655561519	3.792089939
            0.64918309	3.796612978
            0.642804439	3.80178237
            0.63642597	3.807113171
            0.630047256	3.812928438
            0.623668912	3.819067001
            0.617290339	3.825528622
            0.610911952	3.832474947
            0.604533467	3.839905739
            0.598155138	3.84749794
            0.591776797	3.855736494
            0.58539834	3.864298105
            0.579019904	3.873182774
            0.57264159	3.882552147
            0.566262991	3.89208293
            0.559884678	3.902098417
            0.553506117	3.91259861
            0.547127784	3.922936916
            0.540749472	3.933598757
            0.534370921	3.944098711
            0.527992465	3.955406427
            0.521613843	3.966875792
            0.515235278	3.978183508
            0.508856798	3.989329815
            0.502478243	4.000314713
            0.496099768	4.011945248
            0.489721258	4.02357626
            0.48334266	4.035045624
            0.476964007	4.046837807
            0.470585405	4.058792114
            0.464206833	4.070745945
            0.457828255	4.083022594
            0.451449604	4.095299721
            0.445071008	4.107738495
            0.438692341	4.12033844
            0.43231362	4.133100033
            0.425935172	4.145861626
            0.419556692	4.158623219
            0.413178297	4.171546459
            0.406799779	4.184308052
            0.400421399	4.197231293
            0.394042772	4.210154533
            0.387664224	4.222916126
            0.381285829	4.235839367
            0.374907203	4.24860096
            0.368528724	4.261200905
            0.362150406	4.273962498
            0.355772094	4.28656292
            0.34939377	4.299162865
            0.343015136	4.31176281
            0.336636766	4.324201584
            0.330258444	4.336639881
            0.323954472	4.347948074];

%% Interpolations
    U_dis = interp1(data_dis(:,1),data_dis(:,2),x,'linear');
    U_cha = interp1(data_cha(:,1),data_cha(:,2),x,'linear');
    U_avg = (1-chg)*U_dis + chg*U_cha;
end

