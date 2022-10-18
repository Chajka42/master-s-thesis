public class ScoStudentError{
    static float tkoefstd95 = 1.960359441f;
//расчет ошибки по Стьюденту 
    public static void SCO(byte n, double[] charge, int inderx) {

        double midlcrage = 0, errcharge = 0, DRKcharge, tersenterr, DRKchargeSCO, errchargeSCO = 0, tersenterrSCO;
        for (int j = 0; j < inderx-1; j++) {
        charge[j] /= 2002;
        midlcrage += charge[j]; 
        }
        midlcrage /= inderx;
        for (int j = 0; j < inderx-1; j++) {
            errcharge += (charge[j] - midlcrage)*(charge[j] - midlcrage);
            errchargeSCO += (charge[j] - midlcrage)*(charge[j] - midlcrage);
        }
        DRKcharge = (1.0f/(inderx*(inderx-1.0f)));
        DRKchargeSCO = (1.0f/(inderx-1.0f));
        errcharge = Math.sqrt(DRKcharge*errcharge);
        errchargeSCO = Math.sqrt(DRKchargeSCO*errchargeSCO);
        errcharge *= tkoefstd95;
        errcharge /= 2;
        tersenterr = 100/midlcrage*errcharge;
        tersenterrSCO = 100/midlcrage*errchargeSCO;
        /*
        System.out.println("Tok razryada SIN ST  -> " + midlcrageSIN + " +/- " + errchargeSIN + " [мкКл] (" + tersenterrSIN + "%)");
        System.out.println("Tok razryada SIN SCO -> " + midlcrageSIN + " +/- " + errchargeSINSCO + " [мкКл] (" + tersenterrSINSCO + "%)");
        */
        String strokaindex, strokarazmernosti;
              if (n == 1) {strokaindex = "U_ignite"; strokarazmernosti = "kV";} 
        else {if (n == 2) {strokaindex = "Q_charge"; strokarazmernosti = "uC";} 
                     else {strokaindex = "F_neu_pu"; strokarazmernosti = "pu";}}
        System.out.println(strokaindex + " -> " + String.format("%.2f",midlcrage) + " +/- " + String.format("%.2f",errcharge) +
        " ["+strokarazmernosti+"] (" + String.format("%.2f",tersenterr) + "%) // SD -> " + String.format("%.2f",errchargeSCO) +
         " ["+strokarazmernosti+"]; (" + String.format("%.2f",tersenterrSCO) + "%)");
        

    }
}

