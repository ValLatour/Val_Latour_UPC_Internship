package models;

import java.util.*;
public interface Model {
    public String model_name();

    public void dispose();

    public void init_model( double alph,
                            double gam,
                            double delt,
                            String[] F,
                            String[] K,  
                            String[] J,
                            double[][]r,    // r[k][j]
                            List<List<Integer>> Kf, //Kf[f] is the list of indexes of flight f's trajectory in K
                            String[] C,    //array containing all the configurations
                            int [][] Ca, //index of configurations grouped by area control center 
                            int[][] Sl,  // Sl[l] gives an array listing the configurations containing the sector l
                            Integer[] T, //T is the ordered set of time
                            int[][] P, // P = [...;[ind_begin;ind_end];...] periods of T
                            int[][][] Tfkj, // Tfkj[k][j] is an ordered array of time indexes where the flight following k can be in j the indexes are counted from 1
                            int[][]Jfki, // Jfki[k][i] is the ith elementary sector the flight following k passes through 
                            double[] D, // indexed by K
                            double[] E, // indexed by K
                            int [][] TC, // TC[k][j] is the time needed (in number of time steps) for the flight following k to get from j to the next sector
                            int [][] cap,  //cap[l][p] is the capacity of sector l during the period indexed p
                            int [][][] Jktol, // Jktol[l][p][k] represent the first elementary sector for k that function in l in time period p
                            int M);

    public boolean solve(   String log_file,
                            int time_limit_seconds,
                            double heuristics_value,
                            int presolve_value,
                            int[][][] x_res,
                            int[] z_res,
                            int[][] u_res);
    
    public void print_solution( double alph,
                                double gam,
                                double delt,
                                String[] F,
                                String[] K,  
                                String[] J,
                                double[][]r,    // r[k][j]
                                List<List<Integer>> Kf, //Kf[f] is the list of indexes of flight f's trajectory in K
                                String[] C,    //array containing all the configurations
                                int [][] Ca, //index of configurations grouped by area control center 
                                int[][] Sl,  // Sl[l] gives an array listing the configurations containing the sector l
                                Integer[] T, //T is the ordered set of time
                                int[][] P, // P = [...;[ind_begin;ind_end];...] periods of T
                                int[][][] Tfkj, // Tfkj[k][j] is an ordered array of time indexes where the flight following k can be in j the indexes are counted from 1
                                int[][]Jfki, // Jfki[k][i] is the ith elementary sector the flight following k passes through 
                                double[] D, // indexed by K
                                double[] E, // indexed by K
                                int [][] TC, // TC[k][j] is the time needed (in number of time steps) for the flight following k to get from j to the next sector
                                int [][] cap,  //cap[l][p] is the capacity of sector l during the period indexed p
                                int [][][] Jktol, // Jktol[l][p][k] represent the first elementary sector for k that function in l in time period p
                                int M,
                                int[][][] x_res,
                                int[] z_res,
                                int[][] u_res,
                                String output_filename);

}
