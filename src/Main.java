import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.time.format.DateTimeFormatter;  
import java.time.LocalDateTime;    
import java.text.DecimalFormat;

public class Main {
    static DecimalFormat myFormatter = new DecimalFormat("###.##");
    static long start = System.currentTimeMillis();
    int s;
    static int countFiles1;

    //инициализация расчитанной sig(Uuskor) [barn]
    static float [] sigmaparab = new float [2002];

    //t koef student for alpha 0.95 [tkoefstd95.txt] - error value
    static float tkoefstd95 = 1.960359441f;

    static float [] absoluteErrorALL = new float [2002];
    static float [] absoluteErrorSIN = new float [2002];
    static float [] absoluteErrorUNI = new float [2002];
    static float [] absoluteErrorDUO = new float [2002];

    static float [] aErrorTokAll = new float [2002];
    static float [] aErrorTokSIN = new float [2002];
    static float [] aErrorTokUNI = new float [2002];
    static float [] aErrorTokDUO = new float [2002];

    static float [] aNErrorTokSIN = new float [2002];
    static float [] aNErrorTokUNI = new float [2002];
    static float [] aNErrorTokDUO = new float [2002];
    static float [] aNErrorTokALL = new float [2002];

    //нейтронные спектры средние 
    static float [] NEUarraytokmishenSIN = new float [2002];
    static float [] NEUarraytokmishen = new float [2002];
    static float [] NEUarraytokmishenDU = new float [2002];
    static float [] NEUarraytokmishenTOTa = new float [2002];
    
    //относительные нейтронные спектры 
    static float [] otnNEUarraytokmishenSIN = new float [2002];
    static float [] otnNEUarraytokmishen = new float [2002];
    static float [] otnNEUarraytokmishenDU = new float [2002];
    static float [] otnNEUarraytokmishenTOTa = new float [2002];

    //униполярные
    static float [] arrayvolttofile = new float [2002];
    static float [] arraytokmishen = new float [2002];

    //синусоидальные
    static float [] arrayvolttofileSIN = new float [2002];
    static float [] arraytokmishenSIN = new float [2002];

    //двухполярные
    static float [] arrayvolttofileDU = new float [2002];
    static float [] arraytokmishenDU = new float [2002];

    //for total neu
    static float [] arraytokmishenTOTa = new float [2002];

    static int sttype = 0;
    static float sttype2 = 0;
    static float sttype3 = 0;
    static int sttypeUNI = 0;
    static int inderx = 0, indery = 0, inderz = 0, IHJ = 0;
    static float sinneu = 0f, unineu = 0f, duoneu = 0f, sinnn4 = 0f, sinnn8 = 0f, uninn4 = 0f, uninn8 = 0f,
     Gint5 = 0f, Gint7 = 0f, Gint51 = 0f, Gint71 = 0f, GintotSIN = 0f, GintotUNI = 0f, GintotDUO = 0f;
    public static String Ahoii;

    public static void main(String[] args) throws FileNotFoundException, IOException, InterruptedException {

String pathtoparab = AbsPath.PathZero(Ahoii) + "//sigmaparab.txt";
try {
    File file11 = new File("/" + pathtoparab);
    FileReader fr11 = new FileReader(file11);
    try (
    BufferedReader reader11 = new BufferedReader(fr11)) {
        
        String line11 = reader11.readLine();
        int i = 0;
 
        while (line11 != null) {
            i++;
    
            if (i > -1 && i < 2002) { 
               // System.out.println(line);
                String line222 = line11; //обрезаем строчную переменную, оставляя только вольтаж
                try{
                    Float dino11 = Float.parseFloat(line222); 
                    
                    
                    sigmaparab[i] = dino11; 
                   // System.out.println(dino); 
                   
                }
                catch (NumberFormatException ex){
                    ex.printStackTrace(); 
                }
            }
            line11 = reader11.readLine();
        }
    }
} catch (FileNotFoundException e) {
    e.printStackTrace();
} catch (IOException e) {
    e.printStackTrace();
}
        int n = 0, ncase = 1, incase, nn, Gjkcheck = 1;
        

        int countFiles = 0;
        File f= new File(AbsPath.PathZero(Ahoii) + "\\data");
        //System.out.println(AbsPath.PathZero(Ahoii) + "\\data");
        File[] files = f.listFiles();
        for (File file : files) {
            if (file.isFile()) {
                countFiles++;  
            }
        }

     final int FGH = countFiles;

     System.out.println("It looks like we have around: " + countFiles + " files (" + countFiles/3 + " impulses)."+ "\n");
     System.out.println("Predicting calculation time: " +  String.format("%.2f",countFiles/14.5/60) + " minutes."+ "\n");

        float [][] arrayvolt = new float [3][2002]; 

                //создадим массив data base тока мишени
                float [][] datatokm = new float [FGH/3][2002]; //ALL

                //zaryad sbornic integtalov
                double [] chargeSIN = new double [FGH/3];
                double [] chargeUNI = new double [FGH/3];
                double [] chargeDUO = new double [FGH/3];
                double [] chargeALL = new double [FGH/3];

                //tok misheni sbornic integralov
                double [] NEUoSIN = new double [FGH/3];
                double [] NEUoUNI = new double [FGH/3];
                double [] NEUoDUO = new double [FGH/3];
                double [] NEUoALL = new double [FGH/3+1];


                float [] Podzhigmax = new float [FGH/3];

                //сборники типовых импульсов
                float [][] dataSIN = new float [10000][2002];
                float [][] dataUNI = new float [10000][2002];
                float [][] dataDUO = new float [10000][2002];

                float [][] dataTokAll = new float [FGH/3][2002];
                float [][] dataTokSIN = new float [10000][2002];
                float [][] dataTokUNI = new float [10000][2002];
                float [][] dataTokDUO = new float [10000][2002];
                float [] ALLtokrazryda = new float [2002];

                double [][] podgigtopvalue = new double [4][FGH/3];
                
                 float [] VOLTary = new float [FGH/3];
                 float [] VOLTary2 = new float [FGH/3];
                 float [] ticktary = new float [FGH/3];

                    //percentage of sin variation pulse to pulse
                float[] persinvar = new float [FGH/3];
        while (ncase <= (countFiles / 3)){ //главный цикл чтения и обработки
            //генерируем путь к файлам перед циклом чтения
            incase = 0;
            while (incase < 3) {
             String pathto = AbsPath.PathZero(Ahoii) + "\\data\\C0Trace00000.csv"; //адресс папки + название файлов с динамически меняющимся номером

         String strokafirst;  
         //задаем путь
         n = incase + 1; nn = ncase - 1; 
     if (nn < 10000 && nn > 999) {strokafirst = "0" + nn;} 
     else {   if (nn < 1000 && nn > 99) {strokafirst = "00" + nn;} 
     else { if (nn < 100 && nn > 9) { strokafirst = "000" + nn;} 
     else { strokafirst = "0000" + nn;}} }
     pathto = AbsPath.PathZero(Ahoii) + "\\data\\C" + n + "Trace" + strokafirst + ".csv"; 
           try {
            File file = new File("/" + pathto);
            FileReader fr = new FileReader(file);
            try (
            BufferedReader reader = new BufferedReader(fr)) {
                String line = reader.readLine();
                int i = 0;
          //   System.out.println(line);
                while (line != null) {
                    i++;
                   
                    if (i > 5 && i < 2008) { 
                       // System.out.println(line);
                        int index1 = line.indexOf(',');
                        String line2 = line.substring(index1 + 1,line.length()); 
                        try{
                            Float dino = Float.parseFloat(line2); 
                            
                            
                            arrayvolt [incase][i-6] = dino; 
                           // System.out.println(dino); 
                           
                        }
                        catch (NumberFormatException ex){
                            ex.printStackTrace(); 
                        }
                    }
                    line = reader.readLine();
                }
            }
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
             
              incase += 1;
             } //заканчиваем перепись кейса C1-3 0000N
  //обработка C1 
    int Ytick = 0;
    float RR16 = 0, RR6 = 0;
    for (int j =0; j < 1000; j++) {
       if ( RR6 < Math.abs(arrayvolt[0][j]) ) {RR6 = Math.abs(arrayvolt[0][j]);Ytick = j;}
         }  //нашли максимальное значение RR6 по y и Ytick по x
         RR16 = RR6;
         int Startchet6 = Ytick, tfront = 0, Finishtchet6 = Ytick;
         while (RR16 > 0) {RR16 = arrayvolt[0][Finishtchet6+1]; Finishtchet6++;}
         RR16 = RR6;
         while (RR16 > 0) {RR16 = arrayvolt[0][Startchet6-1]; Startchet6--;}
            tfront = Ytick - Startchet6;
    //обработка C2
    float tic1 = 0;
            int tic3 = 0, ticJ2 = 0;
            for (int j = 0; j < arrayvolt[1].length; j++) {
                //пробегает каждое значение внутри импульса
                if (tic1 < arrayvolt[1][j]) {tic1 = arrayvolt[1][j];} //максимальное значение
                if (arrayvolt[1][j] > 0) {ticJ2++;
                } else {if (ticJ2 > 150) {ticJ2 = 0; tic3++;} else {ticJ2 = 0;} } //тип импульса
            }
           //считаем сколько пиков отсекается
            int Startchet = 0;
            float RR1 = 0, RR = 0;
            for (int j = 0; j < 50; j++) {
            RR += Math.abs(arrayvolt[1][j]);  }
                RR1 = RR/50;
                int j = 0;
                while (arrayvolt[1][j] < RR1*30 ) {Startchet = j; j++;}
                j = Startchet;
                int Finishtchet = 0;
                while (arrayvolt[1][j] > 0) {Finishtchet = j; j++;}

                //"сглаживаем"
                Float Nsglag = 75f;
                float [] arrayGLAG = new float [Finishtchet-Startchet + 75];
                arrayGLAG[0] = (float) 0;

                //ищем статистический пик и его границы

                for (int SS = Startchet; SS < Finishtchet + Nsglag; SS++) { 
                    float HESUM = 0f;
                    for (int HE = 0; HE < Nsglag; HE++) {HESUM += arrayvolt[1][SS+HE];}
                    arrayGLAG[SS-Startchet] = HESUM/Nsglag; }
                    
                    float tic4 = 0;
                    int tic4time = 0;
                    for (int IP = 0; IP < Finishtchet - Startchet; IP++) {
                        if (tic4 < arrayGLAG[IP]) {tic4 = arrayGLAG[IP]; tic4time = IP + 37 + Startchet;} }

                        float Supersum = 0f;
                   //     System.out.println(tic4time); 
                   //     System.out.println(Startchet);
                        float [] SEDOT = new float [tic4time - Startchet];
                        for (int IFG = Startchet; IFG < tic4time; IFG++) {
                            SEDOT[IFG - Startchet] = Math.abs(arrayGLAG[IFG - Startchet] - arrayvolt[1][IFG + 37]); 
                        Supersum += SEDOT[IFG - Startchet];
                        }
                        
                        Supersum = Supersum/(tic4time - Startchet); //среднее отклонение 
                        float SEDOTmax = 0; //помогает найти максимальное отклонение 

                        int SDTime = 0;  //время аномалии
                        for (int IFG = Startchet; IFG < tic4time; IFG++) {
                           //если отклонение больше среднего в 3 раза
                           if (SEDOT[IFG - Startchet] > Supersum*3) { 
                               //найти максимальное отклонение и узнать ее время
                            if (SEDOTmax < SEDOT[IFG - Startchet]) {SEDOTmax = SEDOT[IFG - Startchet]; SDTime = IFG;}}

                        }

                        float integral1 = 0f, integral2 = 0f, integral3 = 0f, neutrons = 0f, nn4 = 0f, nn8 = 0f;

                        for (j = 0; j < arrayvolt[1].length; j++) {
                            integral1 += Math.abs(arrayvolt[1][j]);
}
integral1 *= 200;
integral1 /= 100000000;
for (j = 0; j < arrayGLAG.length; j++) {
    integral2 += Math.abs(arrayGLAG[j]);
}
integral2 *= 200;
integral2 /= 10000000;
float integral5 = 0f,  integral7 = 0f, inttot = 0f; 
for (j = 0; j < 2000; j++) {
 //интеграл под ускоряющим ток разряда
 //интеграл под сечением 80keV
 inttot += Math.abs(arrayvolt[1][j]);
 if (j < Startchet6+1056) {if (j > Startchet6) {integral5 += Math.abs(arrayvolt[1][j]);}  else  {integral7+= Math.abs(arrayvolt[1][j]);} }
if (j > Startchet6+1056) {integral7+= Math.abs(arrayvolt[1][j]);}

}
  //обработка C3
  for (j = 500; j < 1500; j++) {
    integral3 += arrayvolt[2][j];
}
integral3 /= 56;
//восстанавливаем спектры в амперы (витки трансформатора и сопротивление шунта)
for (int jen = 0; jen < 2002; jen++) {
arrayvolt[2][jen] /= 56;
arrayvolt[1][jen] *= 200;}

VOLTary[ncase-1] = 0f;
VOLTary2[ncase-1] = 0f;
for (int jf = 0; jf < 2002; jf++) {
    VOLTary[ncase-1] += (float) arrayvolt[2][jf]*sigmaparab[jf]/0.314f;
    VOLTary2[ncase-1] += (float) arrayvolt[2][jf];}
    VOLTary[ncase-1] /= 2002f;
    VOLTary2[ncase-1] /= 2002f;
nn4 += neutrons; nn8 += neutrons;
for (int jf = 0; jf < 2002; jf++) {
    if (Podzhigmax[ncase-1] < arrayvolt[0][jf]) {Podzhigmax[ncase-1] = arrayvolt[0][jf];} 
    for (j = 0; j < 2002; j++) {if (arrayvolt[0][j] > podgigtopvalue[3][ncase-1]) {podgigtopvalue[3][ncase-1] = arrayvolt[0][j];}}
}

//общая статистика
for (j = 0; j < 2002; j++) {
    datatokm[ncase-1][j] =  arrayvolt[2][j];
     dataTokAll[ncase-1][j] = arrayvolt[1][j];
    ALLtokrazryda[j] += arrayvolt[1][j];
    chargeALL[ncase-1]  += Math.abs(arrayvolt[1][j]);
    if (arrayvolt[2][j] > 0) {NEUoALL[ncase]  += Math.abs(arrayvolt[2][j]*sigmaparab[j]/0.314f);}}

//запись C1+C2

String OOI, NESTAB; 

if (SEDOTmax != 0.0f) {NESTAB = "unstable in -> "; NESTAB = NESTAB.substring(0,15) + SDTime;} else {NESTAB = "stable";}
if (tic3 >= 3) {OOI = "sinusoid"; sttype += 1; sinneu += neutrons; inderx += 1; sinnn4 += nn4; sinnn8 += nn8; Gint5 += integral5; 
Gint7 += integral7; GintotSIN += inttot;
        for (j = 0; j < arraytokmishenSIN.length; j++) { arraytokmishenSIN[j] += arrayvolt[2][j];}
                for (j = 0; j < 2002; j++) { arrayvolttofileSIN[j] += arrayvolt[1][j];
                    dataSIN[inderx-1][j] = arrayvolt[2][j];
                     dataTokSIN[inderx-1][j] = arrayvolt[1][j]; 
                    chargeSIN[inderx-1] += Math.abs(arrayvolt[1][j]);
                    if (arrayvolt[2][j] >= 0) {NEUoSIN[inderx-1] += arrayvolt[2][j]*sigmaparab[j]/0.314f; };
                if (arrayvolt[0][j] > podgigtopvalue[0][inderx-1]) {podgigtopvalue[0][inderx-1] = arrayvolt[0][j];}}
                persinvar[inderx-1] = (float) 100/(inderx+indery+inderz)*inderx;

} //тут заканчивается отбор значений (проверка и построение средних)
 else 
{if (tic3 == 2) {OOI = "two peak anomaly"; duoneu += neutrons; inderz += 1; IHJ++; GintotDUO += inttot;
for (j = 0; j < 2002; j++) {arrayvolttofileDU[j] += arrayvolt[1][j];}
        for (j = 0; j < arraytokmishenDU.length; j++) {arraytokmishenDU[j] += arrayvolt[2][j];}
                for (j = 0; j < 2002; j++) { dataDUO[inderz-1][j] = arrayvolt[2][j]; dataTokDUO[inderz-1][j] = arrayvolt[1][j]; 
                    chargeDUO[inderz-1] += Math.abs(arrayvolt[1][j]); 
                    if (arrayvolt[2][j] > 0) {NEUoDUO[inderz-1] += Math.abs(arrayvolt[2][j]*sigmaparab[j]/0.314f);}
                    if (arrayvolt[0][j] > podgigtopvalue[2][inderz-1]) {podgigtopvalue[2][inderz-1] = arrayvolt[0][j];}}} 
else {OOI = "unipolar"; sttypeUNI += 1; unineu += neutrons; indery += 1; uninn4 += nn4; uninn8 += nn8; Gint51 += integral5; 
 Gint71 += integral7; GintotUNI += inttot;
 for (j = 0; j < 2002; j++) {arrayvolttofile[j] += arrayvolt[1][j];}
        for (j = 0; j < arraytokmishen.length; j++) { arraytokmishen[j] += arrayvolt[2][j];}
                for (j = 0; j < 2002; j++) { dataUNI[indery-1][j] = arrayvolt[2][j]; 
                dataTokUNI[indery-1][j] = arrayvolt[1][j]; 
                chargeUNI[indery-1] += Math.abs(arrayvolt[1][j]);
                if (arrayvolt[2][j] > 0) {NEUoUNI[indery-1] += Math.abs(arrayvolt[2][j]*sigmaparab[j]/0.314f);}
                if (arrayvolt[0][j] > podgigtopvalue[1][indery-1]) {podgigtopvalue[1][indery-1] = arrayvolt[0][j];}}}}

//line with information on each impulse
String stroka1 = "Imulse case №" + ncase + ";" + "\n" + "U podzhig: " + "Peak U -> " + RR6 + "; Time peak -> " + Ytick +" ; Start pulse -> " 
+ Startchet6 + "; End pulse -> " + Finishtchet6 + "; Time of front = " + tfront + " ns" +"\n" + 
"I of charge:" + " Type -> " + OOI + "; Peak's Over [1V] detected (type marker) ->" + tic3  
+ "; Umax = " + tic1 + " [V]; Peak U -> " + tic4 + "; Time peak -> " + tic4time +" ; Start pulse -> " 
+ Startchet + "; End pulse -> " + Finishtchet + "; Pulse " + NESTAB + "; INTreal = " + integral1 + "A/c; INTsmooth = " + integral2 + " A/c; " +"\n" +
"U target: INTall = " + integral3 + "A/ms; Neutrons = " + neutrons + "[I_m * sig(U_uskor)];" +"\n";  
  //writing to file !output
 try(FileOutputStream fos = new FileOutputStream(AbsPath.PathZero(Ahoii) + "\\!output.txt", true))
  {
      // перевод строки в байты
      byte[] buffer = stroka1.getBytes();
        
      fos.write(buffer, 0, buffer.length);
  }

int inderx0 = 0, indery0 = 0, inderz0 = 0, tzer = 0;
float uipo;
if (ncase/Gjkcheck == 10) {Gjkcheck++; uipo = (float) 100/(inderx-inderx0+indery-indery0+inderz-inderz0)*(inderx-inderx0); 
   ticktary[tzer] = uipo;
   System.out.println("completed " + ncase + " pulse. (SPHP = " + String.format("%.2f",uipo) + "%)");
   long finish = System.currentTimeMillis();
long elapsed = finish - start;
if (elapsed/1000 < 60) {System.out.println("processing: " + elapsed/1000 + " sec now.");} else 
{System.out.println("processing: " + elapsed/60000 + " min "+ elapsed/1000%60 + " sec now.");}
    
    inderx0 = inderx; indery0 = indery; inderz0 = inderz; tzer += 1;}


ncase += 1;
        }

float p1 = 0f, p2 = 0f, p3 = 0f;
sttype2 = NegativeNeu.PiP(FGH, sttype);
sttype3 = NegativeNeu.PiP(FGH, sttypeUNI);
p1 = NegativeNeu.PiP(FGH, inderx);
p2 = NegativeNeu.PiP(FGH, indery);
p3 = NegativeNeu.PiP(FGH, inderz);

       System.out.println("\n" + " < Oscillatory impulse type - " + String.format("%.2f",p1) + "% > ");
       
       ScoStudentError.SCO((byte) 1, podgigtopvalue[0], inderx);
       ScoStudentError.SCO((byte) 2, chargeSIN, inderx);
       ScoStudentError.SCO((byte) 3, NEUoSIN, inderx);

       System.out.println("\n" + " < Unipolary impulse type - " + String.format("%.2f",p2) + "% > ");

       ScoStudentError.SCO((byte) 1, podgigtopvalue[1], indery);
       ScoStudentError.SCO((byte) 2, chargeUNI, indery);
       ScoStudentError.SCO((byte) 3, NEUoUNI, indery);

       System.out.println("\n" + " < Bipolary impulse type - " + String.format("%.2f",p3) + "% > ");

       ScoStudentError.SCO((byte) 1, podgigtopvalue[2], inderz);
       ScoStudentError.SCO((byte) 2, chargeDUO, inderz);
       ScoStudentError.SCO((byte) 3, NEUoDUO, inderz);

       System.out.println("\n" + " < Impulses of all types - 100% > ");

       ScoStudentError.SCO((byte) 1, podgigtopvalue[3], countFiles/3);
       ScoStudentError.SCO((byte) 2, chargeALL, countFiles/3);
       ScoStudentError.SCO((byte) 3, NEUoALL, countFiles/3);

    //Для интегральных величин в конец !output.txt
        sinneu /= inderx; sinnn4 /= inderx; sinnn8 /= inderx; Gint5 /= inderx; Gint7 /= inderx; GintotSIN /= inderx;
        Gint5 = (Gint5)/10; Gint7 = (Gint7)/10; GintotSIN = (GintotSIN)/10; unineu /= indery; uninn4 /= indery; uninn8 /= indery;
        Gint51 /= indery; Gint71 /= indery; GintotUNI /= indery; Gint51 = (Gint51)/10; Gint71 = (Gint71)/10; GintotUNI = (GintotUNI)/10;
        duoneu /= inderz; GintotDUO /= inderz; GintotDUO /= 10;

        float otnospoterii = 0f, otnosduo = 0f;
        otnospoterii = 1/sinnn4*uninn4;
        otnosduo = 1/sinnn4*duoneu;
           for (int jen = 0; jen < 2002; jen++) {
            arraytokmishenTOTa[jen] += arraytokmishen[jen] + arraytokmishenSIN[jen] + arraytokmishenDU[jen];
           }

for (int jenero = 0; jenero < 2002; jenero++) {
    arrayvolttofile[jenero] /= indery;
     arrayvolttofileSIN[jenero] /= inderx;
      arrayvolttofileDU[jenero] /= inderz;
       ALLtokrazryda[jenero] /= (FGH/3);
       arraytokmishen[jenero] /= indery;
       arraytokmishenSIN[jenero] /= inderx;
       arraytokmishenDU[jenero] /= inderz;
       arraytokmishenTOTa[jenero] /= (inderx+indery+inderz);
   }

//расчет погрешностей для тока на мишени absoluteError
//SIN 
absoluteErrorSIN = StErrCalc.Erro(inderx, dataSIN, arraytokmishenSIN);
aNErrorTokSIN = StErrCalc.ErroNeu(inderx, dataSIN, arraytokmishenSIN, sigmaparab);
//SIN tok
aErrorTokSIN = StErrCalc.Erro(inderx, dataTokSIN, arrayvolttofileSIN);
//UNI
absoluteErrorUNI = StErrCalc.Erro(indery, dataUNI, arraytokmishen);
aNErrorTokUNI = StErrCalc.ErroNeu(indery, dataUNI, arraytokmishen, sigmaparab);
//UNI tok
aErrorTokUNI = StErrCalc.Erro(indery, dataTokUNI, arrayvolttofile);
//DUO
absoluteErrorDUO = StErrCalc.Erro(inderz, dataDUO, arraytokmishenDU);
aNErrorTokDUO = StErrCalc.ErroNeu(inderz, dataDUO, arraytokmishenDU, sigmaparab);
//DUO tok
aErrorTokDUO = StErrCalc.Erro(inderz, dataTokDUO, arrayvolttofileDU);
//ALL 
absoluteErrorALL = StErrCalc.Erro(FGH/3, datatokm, arraytokmishenTOTa);
aNErrorTokALL = StErrCalc.ErroNeu(FGH/3, datatokm, arraytokmishenTOTa, sigmaparab);
//ALL tok
aErrorTokAll = StErrCalc.Erro(FGH/3, dataTokAll, ALLtokrazryda);

//тут на самом деле начинаются нейтроны 24.03
NEUarraytokmishenSIN = NegativeNeu.NegN(arraytokmishenSIN, sigmaparab);
NEUarraytokmishen = NegativeNeu.NegN(arraytokmishen, sigmaparab);
NEUarraytokmishenDU = NegativeNeu.NegN(arraytokmishenDU, sigmaparab);
NEUarraytokmishenTOTa = NegativeNeu.NegN(arraytokmishenTOTa, sigmaparab);

// коэф для относительного выхода 1.7955564f 
//расчет относительных потоков
for (int jen = 0; jen < 2002; jen++) {
    otnNEUarraytokmishenSIN[jen] =  NEUarraytokmishenSIN[jen]/1.7955564f;
     otnNEUarraytokmishen[jen] =  NEUarraytokmishen[jen]/1.7955564f;
      otnNEUarraytokmishenDU[jen] =  NEUarraytokmishenDU[jen]/1.7955564f;
       otnNEUarraytokmishenTOTa[jen] =  NEUarraytokmishenTOTa[jen]/1.7955564f;}

   //запись файлов для графиков 
   System.out.println("\n" + " < Creating Graphs > in directory: " + AbsPath.PathZero(Ahoii) + "\\retail");
   int i1 = 0;
   String AllyNamed = GraphsSaver.Numeric(i1) + "1. U_accel. [kV] Bell-shaped pulse", 
   a1 = "Oscil.", b1 = "Unipol.", c1 = "Bipol.", d1 = "all pulse", e1 = "pulse to pulse";
   GraphsSaver.SaveG(sigmaparab, AllyNamed);
   AllyNamed = "2. Ion I [A] " + a1 ;
   GraphsSaver.SaveG(arraytokmishenSIN, AllyNamed);
   AllyNamed = "3. Ion I [A] " + b1 ;
   GraphsSaver.SaveG(arraytokmishen, AllyNamed);
   AllyNamed = "4. Ion I [A] " + c1 ;
   GraphsSaver.SaveG(arraytokmishenDU, AllyNamed);
   AllyNamed = "5. Ion I [A] " + d1 ;
   GraphsSaver.SaveG(arraytokmishenTOTa, AllyNamed);
   AllyNamed = "6. F_neu [(n/s)/N_tr] " + a1 ;
   GraphsSaver.SaveG(NEUarraytokmishenSIN, AllyNamed);
   AllyNamed = "7. F_neu [(n/s)/N_tr] " + b1 ;
   GraphsSaver.SaveG(NEUarraytokmishen, AllyNamed);
   AllyNamed = "8. F_neu [(n/s)/N_tr] " + c1 ;
   GraphsSaver.SaveG(NEUarraytokmishenDU, AllyNamed);
   AllyNamed = "9. F_neu [(n/s)/N_tr] " + d1 ;
   GraphsSaver.SaveG(NEUarraytokmishenTOTa, AllyNamed);
   AllyNamed = "10. F_neu [pu] " + a1 ;
   GraphsSaver.SaveG(otnNEUarraytokmishenSIN, AllyNamed);
   AllyNamed = "11. F_neu [pu] " + b1 ;
   GraphsSaver.SaveG(otnNEUarraytokmishen, AllyNamed);
   AllyNamed = "12. F_neu [pu] " + c1 ;
   GraphsSaver.SaveG(otnNEUarraytokmishenDU, AllyNamed);
   AllyNamed = "13. F_neu [pu] " + d1 ;
   GraphsSaver.SaveG(otnNEUarraytokmishenTOTa, AllyNamed);
   AllyNamed = "14. F_neu [pu] St.err. " + a1 ;
   GraphsSaver.SaveG(absoluteErrorSIN, AllyNamed);
   AllyNamed = "15. F_neu [pu] St.err. " + b1 ;
   GraphsSaver.SaveG(absoluteErrorUNI, AllyNamed);
   AllyNamed = "16. F_neu [pu] St.err. " + c1 ;
   GraphsSaver.SaveG(absoluteErrorDUO, AllyNamed);
   AllyNamed = "17. F_neu [pu] St.err. " + d1 ;
   GraphsSaver.SaveG(absoluteErrorALL, AllyNamed);
   AllyNamed = "18. I_discharge [A] " + a1 ;
   GraphsSaver.SaveG(arrayvolttofileSIN, AllyNamed);
   AllyNamed = "19. I_discharge [A] " + b1 ;
   GraphsSaver.SaveG(arrayvolttofile, AllyNamed);
   AllyNamed = "20. I_discharge [A] " + c1 ;
   GraphsSaver.SaveG(arrayvolttofileDU, AllyNamed);
   AllyNamed = "21. I_discharge [A] " + d1 ;
   GraphsSaver.SaveG(ALLtokrazryda, AllyNamed);
   AllyNamed = "22. SD I_discharge [A] " + a1 ;
   GraphsSaver.SaveG(aErrorTokSIN, AllyNamed);
   AllyNamed = "23. SD I_discharge [A] " + b1 ;
   GraphsSaver.SaveG(aErrorTokUNI, AllyNamed);
   AllyNamed = "24. SD I_discharge [A] " + c1 ;
   GraphsSaver.SaveG(aErrorTokDUO, AllyNamed);
   AllyNamed = "25. SD I_discharge [A] " + d1 ;
   GraphsSaver.SaveG(aErrorTokAll, AllyNamed);
   AllyNamed = "26. SD F_neu [pu] " + a1 ;
   GraphsSaver.SaveG(aNErrorTokSIN, AllyNamed);
   AllyNamed = "27. SD F_neu [pu] " + b1 ;
   GraphsSaver.SaveG(aNErrorTokUNI, AllyNamed);
   AllyNamed = "28. SD F_neu [pu] " + c1 ;
   GraphsSaver.SaveG(aNErrorTokDUO, AllyNamed);
   AllyNamed = "29. SD F_neu [pu] " + d1 ;
   GraphsSaver.SaveG(aNErrorTokALL, AllyNamed);
   AllyNamed = "30. Oscil. component [%] " + e1 ;
   GraphsSaver.SaveG(persinvar, AllyNamed);
   AllyNamed = "31. F_neu [pu] all type " + e1 ;
   GraphsSaver.SaveG(VOLTary, AllyNamed);
   AllyNamed = "32. Oscil. component [%] per 100 pulse";
   GraphsSaver.SaveG(ticktary, AllyNamed);
   AllyNamed = "33. U_ignite_max [kV] " + d1 ;
   GraphsSaver.SaveG(Podzhigmax, AllyNamed);
   AllyNamed = "34. Integral Ion I [A] all type " + e1 ;
   GraphsSaver.SaveG(VOLTary2, AllyNamed);
      //AllyNamed = "35. Average F_neu [pu] Oscil. " + e1 ;
  // GraphsSaver.SaveG(experneufrompulse, AllyNamed);
     //AllyNamed = "36. Integral Ion I [A] Oscil. " + e1 ;
   //GraphsSaver.SaveG(experneufrompulse2, AllyNamed);

System.out.println(" ");
System.out.println("Calculations completed.");


float YUIGH = 0f, tanukj = 0f;
YUIGH = 1/GintotSIN*GintotDUO;

DateTimeFormatter dtf = DateTimeFormatter.ofPattern("yyyy/MM/dd HH:mm:ss");  
LocalDateTime now = LocalDateTime.now();  
float YFGH = 100 - (sttype2+sttype3);
tanukj = sttype2/100*sinneu + sttype3/100*unineu + YFGH/100*duoneu;
        String strokaend = 
        "\n" + "Total stats:" + "\n" +
         "Sinusoid impulse -> " + sttype2  + "%; "  + "Neutrons_koeff = " + sinneu + " [I_m*sig(U_usk)];" + "\n" +
         "Neu_under_accelerate = " + sinnn4 + "[n/s]; " + "Neu_under_CS_[80KeV] = " + sinnn8 +  "[n/s]; " + "\n" +
         "Ic_total = " + GintotSIN + "[мкКл]; " + "It_under_accelerate = " + Gint5 + "[мкКл]; " + "It_LOST = " + Gint7 + "[мкКл]; " + "\n" + "\n" +


          "Unipolar impulse -> " + sttype3 + "%; " +  "Neutrons_koeff = " + unineu + " [I_m*sig(U_usk)];" + "\n" +
          "Neu_under_accelerate = " + uninn4 + "[n/s]; " + "Neu_under_CS_[80KeV] = " + uninn8 +  "[n/s]; " + "\n" +
          "Ic_total = " + GintotUNI + "[мкКл]; " + "It_under_accelerate = " + Gint51 + "[мкКл]; " + "It_LOST = " + Gint71 + "[мкКл]; " + "\n" + "\n" +


          "Duopolar impulse -> " + YFGH + "%; " +  "Neutrons_[SIN = 1n; DUOpolar = " + otnosduo + " n];" + "\n" + 
          "Ic_total = " + GintotDUO + "[мкКл]; " + "Otnositel'nii tok Sin 1C, Du0 -> "+ YUIGH + " C" + "\n" + "\n" +

          "OTNOSITEL'NII VYHOD (under accel.) -> [SIN = 1n; UNI = " + otnospoterii + "n] " + "\n" + "\n" +
          "TOTAL neu koef = " + tanukj + " [I_m*sig(U_usk)];" + "\n" + 
          "Processing completed " + dtf.format(now) + "\n" +"\n";  
        //запись в файл !output
        try(FileOutputStream fos = new FileOutputStream("C://!output.txt", true))
        {
           
            byte[] buffer = strokaend.getBytes();
              
            fos.write(buffer, 0, buffer.length);
        }
    }
}


