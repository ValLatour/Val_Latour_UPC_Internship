package data;

import java.util.List;

public class Data {
    double alph;
    double gam;
    double delt;
    String[] F;
    String[] K;  
    String[] J;
    double[][]r;    // r[k][j]
    List<Integer>[] Kf; //Kf[f] is the list of indexes of flight f's trajectory in K
    String[] C;    //array containing all the configurations
    int [][] Ca; //index of configurations grouped by area control center 
    int[][] Sl;  // Sl[l] gives an array listing the configurations containing the sector l
    double[] T; //T is the ordered set of time
    int[][] P; // P = [...;[ind_begin;ind_end];...] periods of T
    int[][][] Tfkj; // Tfkj[k][j] is an ordered array of time indexes where the flight following k can be in j the indexes are counted from 1
    int[][]Jfki; // Jfki[k][i] is the ith elementary sector the flight following k passes through 
    double[] D; // indexed by K
    double[] E; // indexed by K
    int [][][] TC; // TC[k][j][jp] is the time needed (in number of time steps) for the flight following k to get from j to jp (sectors)
                    // if k does not go from j to jp let TC[k][j][jp] be superior to T.length
    int [][] cap;  //cap[l][p] is the capacity of sector l during the period indexed p
    int [][][] Jktol; // Jktol[l][p][k] represent the first elementary sector for k that function in l in time period p
    int M;

    public double getAlph() {
        return alph;
    }
    public void setAlph(double alph) {
        this.alph = alph;
    }
    public String[] getC() {
        return C;
    }
    public void setC(String[] c) {
        C = c;
    }
}
