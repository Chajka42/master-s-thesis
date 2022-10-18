public class NegativeNeu {

    public static float[] NegN(float[] arraytokmishenSIN, float[] sigmaparab) {
        float [] NEUarraytokmishenSIN = new float [2002];
        
        for (int jen = 0; jen < 2002; jen++) {
            if (arraytokmishenSIN[jen] <= 0){NEUarraytokmishenSIN[jen] = 0f;} else {
            NEUarraytokmishenSIN[jen] = arraytokmishenSIN[jen] * sigmaparab[jen];}
           }
        return NEUarraytokmishenSIN;
 
}
    public static float PiP(int Counts, int stt) {
        float pp = 0f;
        pp = 100/(Counts/3)*stt;
        return pp;
    }

}