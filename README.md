# COVID-19 Salgin Modellemesi

## Kompartmanlar

| Notasyon   | Tanim | 
| ------------- |-------------| 
| S | (Sensitive) Henuz enfekte olmamis, hastaliga acik bireyler. | 
| E | (Exposed) Enfekte olmus, fakat henuz bulasici olmayan bireyler. | 
| I | (Infected) Populasyondaki bulasici bireyler. | 
| H | (Hospitilized) Hastanedeki (servis) bireyler. | 
| ICU | (Intensive Care Unit) Yogun bakimdaki bireyler. | 
| R | (Recovered) Iyilesmis ve bagisiklik kazanmis toplam birey sayisi. | 
| X | Vefat etmis toplam birey sayisi. | 
| C | Toplam vaka sayisi (semptomatik + semptomatik) | 
| C<sub>1</sub> | Toplam vaka sayisi (eve gonderilmis, semptomatik) | 
| C<sub>2</sub> | Toplam vaka sayisi (eve gonderilmis,asemptomatik) | 
| R<sub>c</sub> | Kaydedidilen toplam iyilesen sayisi | 

## Parametreler

| Notasyon   | Tanim | 
| ------------- |-------------| 
| R<sub>0</sub> | Temel bulastirma katsayisi | 
| &beta; | Temas orani | 
| r<sub>L</sub> | Karantina etkisi | 
| m<sub>L</sub> | Karantinanin etkisini gosterme egimi| 
| s<sub>L</sub> | Karantinanin etkisini gosterme gecikmesi | 
| 1/τ | Inkubasyon suresi | 
| 1/γ<sub>s</sub> | Bulastiricilik suresi | 
| 1/γ<sub>H</sub> | Hastanede (servis) kalma suresi | 
| 1/γ<sub>ICU</sub> | ICU da kalma suresi | 
| 1/ε<sub>H</sub> | Tanisi konan semptomatiklerin hasteneye kaldirilma orani | 
| 1/ε<sub>H2I</sub> | Sevisten ICU'ya transfer orani | 
| 1/ε<sub>H2x</sub> | Servisten vefat orani | 
| 1/ε<sub>I2x</sub> | ICU'dan vefat orani | 
| r<sub>d</sub><sup>s</sup>| Semptomatiklerin detekte edilme orani | 
| r<sub>d</sub><sup>a</sup>| Asemptomatiklerin detekte edilme orani | 
| r<sub>r</sub><sup>s</sup>| Detekte edilmis semptomatiklerin iyilesme orani |  
| r<sub>r</sub><sup>a</sup>| Detekte edilmis asemptomatiklerin iyilesme orani |  


## Yardimci Fonksiyonlar

Karantina etkisi (r<sub>lock</sub>(t)): Zamana bagli sigmoidal bir fonskiyon olarak modellenmistir. <br/>
r<sub>lock</sub>(t) = r<sub>L</sub> + (1-r<sub>L</sub>) / (1+ exp[m<sub>L</sub>(t-t<sub>L</sub>-s<sub>L</sub>)])

Gevseme etkisi (r<sub>relax</sub>(t)) : Zamana bagli sigmoidal bir fonskiyon olarak modellenmistir. <br/>
r<sub>relax</sub>(t) = r<sub>R</sub> + 1/(1/(1-r<sub>R</sub>) + exp[-m<sub>R</sub>(t-t<sub>R</sub>-s<sub>R</sub>)])

Test kapasitesi etkisi (r<sub>test</sub>(t)) : Zamana bagli ucuncu dereceden bir polinom olarak modellenmistir. Asemptomatiklerin detekte olma oranini etkiler.


## Diferansiyel Denklemler
---------------------
dS / dt = 

h<sub>&theta;</sub>(x) = &theta;<sub>o</sub> x + &theta;<sub>1</sub>x
