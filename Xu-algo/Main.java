//import java.io.*;                 //should already be imported with DataCreator
import java.util.*;               //should already be imported with DataCreator

import models.*;
import data.CapacitySelector;
import data.DataManager;
import data.Minimum_capacity;

public class Main {
    public static Model model_choice(String name){
        switch (name) {
            case "SC_DCB_I":
                return new SC_DCB_I_Model();
            case "Benchmark":
                return new Benchmark_Model();
            default:
                System.out.println("The model "+name+" has not been implemented. If it has, please check for spelling errors.");
                return new SC_DCB_I_Model();
        }
    }

    public static CapacitySelector selector_choice(String name){
        switch (name) {
            case "Minimum_capacity":
                return new Minimum_capacity();
            default:
                System.out.println("The model "+name+" has not been implemented. If it has, please check for spelling errors.");
                return new Minimum_capacity();
        }
    }

    public static void main(String[] args){
        String clusters_file = "";
        String mcx_file = "";
        String trajectory_cost_file = "";
        String configuration_file = "";
        String airspace_file = "";
        String sector_cost_file = "";
        String output_solution_file = "";
        String log_file = "";
        int time_limit_seconds;
        double heuristics_value;
        int presolve_value;
        Model mod;
        CapacitySelector selector;

        if(args.length>0)
        {
            System.out.println("Reading files from args");
            clusters_file = args[0];
            mcx_file = args[1];
            trajectory_cost_file = args[2];
            configuration_file = args[3];
            airspace_file = args[4];
            sector_cost_file = args[5];
            output_solution_file = args[6];
            log_file = args[7];
            time_limit_seconds = Integer.parseInt(args[8]);
            heuristics_value = Double.parseDouble(args[9]) ;
            presolve_value = Integer.parseInt(args[10]);
            mod = model_choice(args[11]);
            selector = selector_choice(args[12]);

            System.out.println("mcx_file: " + mcx_file);
            System.out.println("trajectory_cost_file: " + trajectory_cost_file);
            System.out.println("configuration_file: " + configuration_file);
            System.out.println("airspace_file: " + airspace_file);
            System.out.println("sector_cost_file: " + sector_cost_file);
            System.out.println("output_solution_file: " + output_solution_file);
            System.out.println("log_file: " + log_file);
            System.out.println("Time limit: " + time_limit_seconds + " s");
            System.out.println("Heuristics value: " + heuristics_value);
            System.out.println("Presolve value: " + presolve_value);
        }
        else
        {
            mod = new SC_DCB_I_Model();
            selector = new Minimum_capacity();
            clusters_file = "/mnt/4b80ee63-1fe6-4fa9-92cd-8846477d3287/Data/Marc_PhD/Complexity_generation/model_III/test_annealing/clusters.txt";
            mcx_file = "/mnt/4b80ee63-1fe6-4fa9-92cd-8846477d3287/Data/Marc_PhD/Complexity_generation/model_III/complexity_files/clean/study2_2016-07-28T00:00_2016-07-28T02:00_1469664000.mcx";
            //mcx_file = "/mnt/4b80ee63-1fe6-4fa9-92cd-8846477d3287/Data/Marc_PhD/Complexity_generation/model_III/complexity_files/clean/study2_2016-07-28T02:00_2016-07-28T04:00_1469664000.mcx";
            trajectory_cost_file = "/mnt/4b80ee63-1fe6-4fa9-92cd-8846477d3287/Data/Marc_PhD/Complexity_generation/model_III/summary_files/trajectory_costs.csv";
            configuration_file = "/mnt/4b80ee63-1fe6-4fa9-92cd-8846477d3287/Data/Marc_PhD/ddr2_inputs/airspace/20160728/ECAC/Configuration_1608.cfg";
            airspace_file = "/mnt/4b80ee63-1fe6-4fa9-92cd-8846477d3287/Data/Marc_PhD/Complexity_generation/model_III/airspace/FABEC.spc";
            sector_cost_file = "/mnt/4b80ee63-1fe6-4fa9-92cd-8846477d3287/Data/Marc_PhD/Complexity_generation/model_III/airspace/cost_100.csv";
            output_solution_file = "/mnt/4b80ee63-1fe6-4fa9-92cd-8846477d3287/Data/Marc_PhD/Complexity_generation/model_III/test_annealing/test_MILP_10.txt";
            log_file = "/mnt/4b80ee63-1fe6-4fa9-92cd-8846477d3287/Data/Marc_PhD/Complexity_generation/model_III/test_annealing/log_test_MILP10.txt";
            time_limit_seconds = Integer.MAX_VALUE;
            heuristics_value = 0.05; //Default 0.05 https://www.gurobi.com/documentation/9.0/refman/heuristics.html
            presolve_value = -1; //Default -1 (automatic) https://www.gurobi.com/documentation/9.5/refman/presolve.html#parameter:Presolve

            //Toy example
//            clusters_file = "/mnt/4b80ee63-1fe6-4fa9-92cd-8846477d3287/Data/Marc_PhD/Complexity_generation/model_III/toy_example_files/clusters.txt";
//            mcx_file = "/mnt/4b80ee63-1fe6-4fa9-92cd-8846477d3287/Data/Marc_PhD/Complexity_generation/model_III/toy_example_files/complexity.mcx";
//            trajectory_cost_file = "/mnt/4b80ee63-1fe6-4fa9-92cd-8846477d3287/Data/Marc_PhD/Complexity_generation/model_III/toy_example_files/cost_table.csv";
//            configuration_file = "/mnt/4b80ee63-1fe6-4fa9-92cd-8846477d3287/Data/Marc_PhD/Complexity_generation/model_III/toy_example_files/configuration.cfg";
//            airspace_file = "/mnt/4b80ee63-1fe6-4fa9-92cd-8846477d3287/Data/Marc_PhD/Complexity_generation/model_III/toy_example_files/airspace.spc";
//            sector_cost_file = "/mnt/4b80ee63-1fe6-4fa9-92cd-8846477d3287/Data/Marc_PhD/Complexity_generation/model_III/toy_example_files/cost_sectors_toy.csv";
//            output_solution_file = "/mnt/4b80ee63-1fe6-4fa9-92cd-8846477d3287/Data/Marc_PhD/Complexity_generation/model_III/toy_example_files/output_test_2.txt";
//            log_file = "/mnt/4b80ee63-1fe6-4fa9-92cd-8846477d3287/Data/Marc_PhD/Complexity_generation/model_III/toy_example_files/log_test_2.txt";
        }

        //Config parameters
        double max_complx = 500.0;
        int period_duration = 15; //[min] 10 15
        // Penalty cost for different configurations 10 0
        double phi = 0; //10 0
        // Cost of opening sector s during the time t - Index: {t,s}
        //double theta_value = 100.0;
        //For linearisation -> Big value
        int M   = 1000;
        //int M_h = 1000;
        double threshold_value = 1300.0; //1300
        //double threshold_value = 2.0; //Toy example

        //data structure
        double alph;            //cost of delay
        double gam;             //cost of fuel
        double delt;            //cost of opening configuration/sector
        String[] F;             //flights
        String[] K;             //trajectories
        String[] J;             // elementary sectors
        double[][]r;            // r[k][j]; estimated time of departure
        List<List<Integer>> Kf;     //Kf[f] is the list of indexes of flight f's trajectory in K
        String[] C;             //array containing all the configurations
        int [][] Ca;            //index of configurations grouped by area control center 
        int[][] Sl;             // Sl[l] gives an array listing the configurations containing the sector l
        Integer[] T;             //T is the ordered set of time (no jump)
        int[][] P;              // P = [...;[ind_begin;ind_end];...] periods of T (with their index of beginning/ending)
        int[][][] Tfkj;         // Tfkj[k][j] is an ordered array of time indexes where the flight following k can be in j the indexes are counted from 1
        int[][]Jfki;            // Jfki[k][i] is the ith elementary sector the flight following k passes through 
        double[] D;             // D[k] consumption of fuel for trajectory k 
        double[] E;             // E[k] cost of selecting trajectory k
        int [][] TC;            // TC[k][j] is the time needed (in number of time steps) for the flight following k to get from j to the next sector
        int [][] cap;           //cap[l][p] is the capacity of sector l during the period indexed p
        int [][][] Jktol;       // Jktol[l][p][k] represent the first elementary sector for k that function in l in time period p

         // Loading the data
        //load_toy_example_data(F,K,Cluster_set,C,S,Se,T,P,times_in_p,ind_times_in_p,Kf,Eps,C_s,conf_in_cluster,sectors_in_conf,delta_cost,threshold,C_kk,B,theta_value,theta_period);

        DataManager.initialise_data(max_complx,period_duration,clusters_file, mcx_file, trajectory_cost_file, configuration_file,airspace_file,sector_cost_file,F,K,Cluster_set,C,S,Se,T,P,Kf,Eps,C_s,conf_in_cluster,cluster_conf,sectors_in_conf,delta_cost,threshold_value,threshold,C_kk,B,theta_period,selector);
        DataManager.print_data_summary(F,K, J, P, C,T,Sl);

        M = (int) (max_complx*F.length);
        //M_h = (int) (max_complx*F.size()*F.size());

        System.out.println("M = " + M);
        //System.out.println("M_h = " + M_h);
        System.out.println();

        //Strutures for storing output
        int[][][] x_res;        //delay allocation
        int[] z_res;            //trajectory selection
        int[][] u_res;          //configuration selection

        mod.init_model(alph, gam, delt, F, K, J, r, Kf, C, Ca, Sl, T, P, Tfkj, Jfki, D, E, TC, cap, Jktol, M);
        mod.solve(log_file, time_limit_seconds, heuristics_value, presolve_value, x_res, z_res, u_res);
        mod.print_solution(alph, gam, delt, F, K, J, r, Kf, C, Ca, Sl, T, P, Tfkj, Jfki, D, E, TC, cap, Jktol, M, x_res, z_res, u_res, log_file);
        mod.dispose();

        System.out.println("Finished");
    }
}
