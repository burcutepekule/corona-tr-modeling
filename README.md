# COVID-19 Salgin Modellemesi

## Kompartmanlar

| Notasyon   | Tanim | 
| ------------- |-------------| 
| S | (Susceptible) Henuz enfekte olmamis, hastaliga acik bireyler. | 
| E | (Exposed) Enfekte olmus, fakat henuz bulasici olmayan bireyler. | 
| I | (Infected) Populasyondaki bulasici bireyler. | 
| H | (Hospitilized) Hastanedeki (servis) bireyler. | 
| ICU | (Intensive Care Unit) Yogun bakimdaki bireyler. | 
| R | (Recovered) Iyilesmis ve bagisiklik kazanmis toplam birey sayisi. | 
| X | Vefat etmis toplam birey sayisi. | 
| C | Toplam vaka sayisi (semptomatik + asemptomatik) | 

## Parametreler

| Notasyon   | Tanim | 
| ------------- |-------------| 
| N | Toplam nufus | 
| R<sub>0</sub> | Temel bulastirma katsayisi | 
| r<sub>L</sub> | Karantina etkisi | 
| m<sub>L</sub> | Karantinanin etkisini gosterme egimi| 
| s<sub>L</sub> | Karantinanin etkisini gosterme gecikmesi | 
| 1/τ | Inkubasyon suresi | 
| 1/γ<sub>s</sub> | Bulastiricilik suresi | 
| 1/γ<sub>H</sub> | Hastanede (servis) kalma suresi | 
| 1/γ<sub>ICU</sub> | ICU da kalma suresi | 
| ε<sub>H2I</sub> | Sevisten ICU'ya transfer orani | 
| ε<sub>H2x</sub> | Servisten vefat orani | 
| ε<sub>I2x</sub> | ICU'dan vefat orani | 
| r<sub>d</sub><sup>s</sup>| Semptomatiklerin detekte edilme orani (hastaneye kaldirildiklari varsayilir) | 
| r<sub>d</sub><sup>a</sup>| Asemptomatiklerin detekte edilme orani (hastaneye kaldirilmadiklari varsayilir)| 


## Yardimci Fonksiyonlar

Karantina etkisi (r<sub>lock</sub>(t)): Zamana bagli sigmoidal bir fonskiyon olarak modellenmistir. <br/>
r<sub>lock</sub>(t) = r<sub>L</sub> + (1-r<sub>L</sub>) / (1+ exp[m<sub>L</sub>(t-t<sub>L</sub>-s<sub>L</sub>)])

Gevseme etkisi (r<sub>relax</sub>(t)) : Zamana bagli sigmoidal bir fonskiyon olarak modellenmistir. <br/>
r<sub>relax</sub>(t) = r<sub>L</sub> + 1/(1/(r<sub>end</sub>-r<sub>L</sub>) + exp[-m<sub>R</sub>(t-t<sub>R</sub>-s<sub>R</sub>)])

Bulastirma carpani (coeff<sub>R</sub>(t)) : Zamana bagli bulastirma katsayisi carpani (karantina ya da gevseme durumuna gore farkli degerler alir)

Karantina durumunda coeff<sub>R</sub>(t) = (1/N)r<sub>lock</sub>(t)<br/>
Gevseme durumunda coeff<sub>R</sub>(t) = (1/N)r<sub>relax</sub>(t)

Test kapasitesi etkisi (r<sub>test</sub>(t)) : Zamana bagli lognormal bir fonksiyon olarak modellenmistir. Asemptomatiklerin detekte olma oranini etkiler. 

#### Test kapasitesinin parametrizasyonu 
<img src="https://github.com/burcutepekule/corona-tr-modeling/blob/master/OUT_17_Jun_2020/FIGS/TESTS_estimate.png" width="625" height="525">

## Diferansiyel Denklemler

dS(t) / dt    = - coeff<sub>R</sub>(t)R<sub>0</sub>γ<sub>s</sub>S(t)I(t)<br/>
dE(t) / dt    = + coeff<sub>R</sub>(t)R<sub>0</sub>γ<sub>s</sub>S(t)I(t) - τE<br/>
dI(t) / dt    = + τE - γ<sub>s</sub>I(t) <br/>
dH(t) / dt    = + r<sub>d</sub><sup>s</sup>γ<sub>s</sub>I(t) - γ<sub>H</sub>H(t) <br/>
dICU(t) / dt  = + γ<sub>H</sub>ε<sub>H2I</sub>H(t) - γ<sub>ICU</sub> ICU(t) <br/>
dR(t) / dt    = + γ<sub>H</sub>(1-ε<sub>H2I</sub>-ε<sub>H2x</sub>)H(t) + γ<sub>ICU</sub>(1-ε<sub>I2x</sub>)ICU(t) + (1-r<sub>d</sub><sup>s</sup>)γ<sub>s</sub>I(t) <br/>
dX(t) / dt    = + γ<sub>H</sub>ε<sub>H2x</sub>H(t) +  γ<sub>ICU</sub>ε<sub>I2x</sub>ICU(t) <br/>
dC(t) / dt    = + (r<sub>d</sub><sup>s</sup> + r<sub>test</sub>(t)r<sub>d</sub><sup>a</sup>)γ<sub>s</sub>I(t)<br/>

## Kompartmanlar arasi gecisi anlatan model diyagrami

<img src="https://github.com/burcutepekule/corona-tr-modeling/blob/master/MODEL_DESCP.png" width="575" height="301">

## Parametre kestirimi icin kullanilan sinyaller

- Gunluk vaka sayisi toplam vaka sayisi kompartmaninin (C(t)) birinci turevine (Gunluk C(t)) oturtulur. <br/>
- Gunluk vefat sayisi toplam vefat sayisi kompartmaninin (X(t)) birinci turevine (Gunluk X(t)) oturtulur. <br/>
- Yogun bakimdaki hasta sayisi ICU kompartmanina oturtulur. <br/>
- Gunluk iyilesen sayisinin kaydinin nasil tutulduguna dair tam bir bilgi(m) olmadigi icin simdilik parametre kestirimine dahil edilmemistir. 

## Veri oturtma 

- Sinyaller uzerindeki hata Negatif Binomial dagilimina gore modellenir, sacilim parametresi her sinyale ozel olarak kestirilir.
- Toplam yerine gunluk veri noktalarinin kullanilmasinin sebebi korele olmayan veri noktalari elde etmektir.
- Verilerin oturtulmasi surecinde Hamilton MCMC algoritmasi kullanilmistir.

## Kod dosyalari
### Gevseme goz onunde bulundurlmadan

- Kosulmasi gereken ana kod dosyasi : ``model_fitting.R``. 
- Verinin on islenmesi icin kosulan kod dosyasi : ``prepare_data.R``. 
- Gerekli fonksiyonlarin cagirildigi kod dosyasi : ``setup.R``. 
- MCMC sonuclarinin analiz edildigi kod dosyasi : ``analysis_chains.R``. 
- Model ciktilarinin gorsellestirildigi kod dosyasi : ``analysis_plots.R``. 
- Gunluk vaka ve seri aralik dagilimi kullanilarak efektif R'in hesaplandigi kod dosyasi : ``analysis_RE.R``. 
- Gunluk test sayisi ve pozitiflik oranina gore duzeltilmis efektif R'in hesaplandigi kod dosyasi : ``analysis_RE_CORR.R``. 
- Normalize edilmis test kapasitesinin lognormal fonskiyona oturtularak gorsellestirildigi kod dosyasi : ``analysis_TESTS.R``. 

### Gevseme goz onunde bulundurularak

- Kosulmasi gereken ana kod dosyasi : ``model_fitting_RELAX_REND.R``. 
- Verinin on islenmesi icin kosulan kod dosyasi : ``prepare_data.R``. 
- Gerekli fonksiyonlarin cagirildigi kod dosyasi : ``setup.R``. 
- MCMC sonuclarinin analiz edildigi kod dosyasi : ``analysis_chains_RELAX_REND.R``. 
- Model ciktilarinin gorsellestirildigi kod dosyasi : ``analysis_plots_RELAX_REND.R``. 
- Gunluk vaka ve seri aralik dagilimi kullanilarak efektif R'in hesaplandigi kod dosyasi : ``analysis_RE.R``. 
- Gunluk test sayisi ve pozitiflik oranina gore duzeltilmis efektif R'in hesaplandigi kod dosyasi : ``analysis_RE_CORR.R``. 
- Normalize edilmis test kapasitesinin lognormal fonskiyona oturtularak gorsellestirildigi kod dosyasi : ``analysis_TESTS.R``.

## Kullanim

- Model ciktilarinin elde edilmesi, gorsellestirilmesi, ve MCMC sonuclarinin analizi icin yalnizca ``model_fitting.R`` dosyasinin kosulmasi yeterlidir. Bu kod dosyasi ``~/DATA`` alt klasorunden ``~/CORONA_TR.csv`` adli Turkiye COVID-19 verilerinin tutuldugu csv dosyasini okuyarak veriyi modele oturtur.

- Model ciktilari kodlarin kosuldugu gunde gore isimlendirilen alt klasore ``~/OUT_<gun>_<ay>_<yil>`` adi altinda kaydedilir. Bu klasorun altinda 3 alt klasor daha olusturulur (``~/CSVS``, ``~/FIGS``, ``~/RDATA``) ve model ciktilari turlerine ait olan klasorlere kaydedilir.  

- Hamiltonian MCMC hesaplama acisindan kaynak kullanimi yuksek olan bir algoritma oldugu icin kosma suresi kullanilan isinma evresi (warmup), iterasyon sayisi (iter), ve zincir (chains) sayisina gore degisecektir. Bu parametreler ``model_fitting.R`` dosyasinin icinden degistirilebilir, ve kosma suresi isinma evresi ve iterasyon sayisi kisaltilarak azaltilabilir. Fakat bu kisaltma sonuclarin guven araligini ve sonsal dagilimlarin yakinsama performasini etkileyebilir.

## Projeksiyonlar
### Son Guncellenme Tarihi : 17 Haziran 2020

### Populasyon Degerleri

Gevseme analizi : Gevseme 1 Haziran merkezli, karantina etkisine benzeyen bir sigmoidal fonksiyon olarak modellenmistir.  r<sub>end</sub> gevsemenin erisecegi maksimum katsayiyi temsil etmektedir. 

![Folder Structure](https://github.com/burcutepekule/corona-tr-modeling/blob/master/OUT_17_Jun_2020/FIGS/figure_Re_mrelax_1.png)

![Folder Structure](https://github.com/burcutepekule/corona-tr-modeling/blob/master/OUT_17_Jun_2020/FIGS/figure_all_mrelax_1.png)

### Efektif R
Efektif R iki sekilde hesaplanabilir. Birincisi model ciktisina gore yapilan efektif R hesaplamasidir, fakat bu hesaplama oturtulan egriye gore yapildigi icin son gunlerdeki vaka sayisindaki dalgalanmalarin etkisini icermemektedir. 

Dalgalanmalari hesaba katan, daha hassas bir efektif R hesaplamasi icin ``analysis_RE.R`` kod dosyasi kullanilabilir. Bu dosya modele gore hesaplanan seri aralik dagilimini ve gunluk vaka sayisindaki degisimi kullanarak zamana bagli bir efektif R kestirimi yapar. 

![Folder Structure](https://github.com/burcutepekule/corona-tr-modeling/blob/master/OUT_17_Jun_2020/FIGS/RE_estimate.png)

#### Referans icin : http://doi.org/10.5281/zenodo.3835635







