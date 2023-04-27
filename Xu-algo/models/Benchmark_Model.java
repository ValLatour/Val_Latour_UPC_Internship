
package models;

import gurobi.*;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import static java.lang.Math.max;
import static java.lang.Math.min;


public class Benchmark_Model implements Model {

    GRBEnv env;
    GRBModel model;
    GRBVar[][][] x;
    GRBVar[] z;
    int nt;
    int nj;
    int nk;

    public String model_name(){return "Benchmark";}

    public Benchmark_Model(){
        //Model
        env = new GRBEnv();
        model = new GRBModel(env);
        model.set(GRB.StringAttr.ModelName, "Benchmark model");
    }

    @Override
    public void init_model(  double alph,
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
                        int [][] TC, // TC[k][j][jp] is the time needed (in number of time steps) for the flight following k to get from j to jp (sectors)
                                        // if k does not go from j to jp let TC[k][j][jp] be superior to T.length
                        int [][] cap,  //cap[l][p] is the capacity of sector l during the period indexed p
                        int [][][] Jktol, // Jktol[l][p][k] represent the first elementary sector for k that function in l in time period p
                        int M)
        {
        int nt = T.length;
        int nj = J.length;
        int nk = K.length;

        try{
            long start_model_init_time = System.currentTimeMillis();

            // Create decision variable for delay selection
            System.out.println("Creating decision variable x");
            x = new GRBVar[nk][nj][nt+1];                       // +1 because we consider the time step 0
            for (int k=0;k<nk;k++){
                for (int j = 0; j <nj; j++){
                    for (int t = 0; t < nt+1; t++){
                        x[k][j][t] = model.addVar(0,1,0,GRB.BINARY,                 //not indexed by f because k already olds this information
                            "x_k"+K[k]+"_j"+J[j]+"_t"+t);
                    }
                }
            }

            //Create decision variable for trajectory selection
            System.out.println("Creating decision variable z");
            z = new GRBVar[nk];
            for (int k = 0 ; k<nk ; k++){
                z[k]=model.addVar(0,1,0,GRB.BINARY,"z_"+K[k]);
            }

            //Create intermediary variable for entry count
            System.out.println("Creating intermediary variable for entry count");
            GRBVar[][][] entry = new GRBVar[nk][nj][nt];
            for (int k=0;k<nk;k++){
                for (int j = 0; j <nj; j++){
                    for (int t = 0; t < nt; t++){
                        entry[k][j][t] = model.addVar(0,1,0,GRB.BINARY,                 //not indexed by f because k already olds this information
                            "entry_k"+K[k]+"_j"+J[j]+"_t"+t);
                    }
                }
            }

            //Objective function
            System.out.println("Creating objective function");
            GRBLinExpr obj = new GRBLinExpr();
            for(int f=0; f<nf;f++){
                for(int k:Kf.get(f)){
                    //cost of delay
                    int j = Jfki[k][0];
                    for (int ti : Tfkj[k][j]){
                        obj.addTerm(alph*(T[ti-1]-r[k][j]),entry[k][j][ti]);
                    }
                    //cost of trajectory
                    obj.addTerm(gam*D[k]+E[k],z[k]);
                }
            }
            model.setObjective(obj,GRB.MINIMIZE);

            //Constraint counting entries
            System.out.println("Creating constraint counting entries");
            for (int k=0;k<nk;k++){
                for (int j = 0; j <nj; j++){
                    for (int t = 1; t < nt+1; t++){
                        GRBLinExpr e_counting_x = new GRBLinExpr();
                        e_counting_x.addTerm(1,x[k][j][t]);
                        e_counting_x.addTerm(-1,x[k][j][t-1]);
                        GRBLinExpr e_counting_entry = new GRBLinExpr();
                        e_counting_entry.addTerm(1,entry[k][j][t]);
                        model.addConstr(e_counting_x,GRB.EQUAL,e_counting_entry,
                            "Count entry_"+K[k]+"_"+J[j]+"_t"+t);
                    }
                }
            }

            //Constraint "exactly 1 trajectory per flight"
            System.out.println("Creating constraint one trajectory per flight");
            for(int f=0; f<nf;f++){
                GRBLinExpr Itraj_f = new GRBLinExpr() ;
                for (int k : Kf.get(f)){
                    Itraj_f.addTerm(1,z[k]);
                }
                model.addConstr(Itraj_f,GRB.EQUAL,1,"Only_one "+F[f]);
            }

            //Constraint minimal & maximal arrival time in sector
            System.out.println("Creating constraints minimal & maximal arrival time in sector");
            for(int f=0; f<nf;f++){
                for (int k:Kf.get(f)){
                    for (int j = 0; j <nj; j++){
                        int before = Tfkj[k][j][0]-1;                   //index of time before possibility of arrival
                        int last = Tfkj[k][j][Tfkj[k][j].length -1];    //last index of time for arrival
                        GRBLinExpr min_arrival = new GRBLinExpr();
                        min_arrival.addTerm(1,x[k][j][before]);
                        model.addConstr(min_arrival,GRB.EQUAL,0,
                            "Minimal_arrival_time k"+K[k]);
                        GRBLinExpr max_arrival = new GRBLinExpr();
                        max_arrival.addTerm(1,x[k][j][last]);
                        model.addConstr(max_arrival,GRB.EQUAL,z[k],
                            "Maximal_arrival_time k"+K[k]);
                    }
                }
            }

            //Constraint "if we previously arrived, it is true in the future"
            System.out.println("Creating constraints definitive arrival");
            for(int f=0; f<nf;f++){
                for (int k: Kf.get(f)){
                    for (int j = 0; j <nj; j++){
                        for (int t = 1; t < nt+1; t++){
                            GRBLinExpr arrived = new GRBLinExpr();
                            arrived.addTerm(1,x[k][j][t]);
                            arrived.addTerm(-1,x[k][j][t-1]);
                            model.addConstr(arrived,GRB.GREATER_EQUAL,0,
                                "Def_arrival "+K[k]+" to "+J[j]);
                        }
                    }
                }
            }

            //Constraint "no airborne delay"
            System.out.println("Creating constraints no airborne delay");
            for(int f=0; f<nf;f++){
                for (int k: Kf.get(f)){
                    int jactu = Jfki[k][0];
                    for (int i=1;i<Jfki[k].length;i++){
                        int jnext = Jfki[k][i];
                        int travel_time = TC[k][jactu];
                        for (int t = 0; t < nt+1 - travel_time; t++){
                            GRBLinExpr arrived = new GRBLinExpr();
                            arrived.addTerm(1,x[k][jnext][t+travel_time]);
                            arrived.addTerm(-1,x[k][jactu][t]);
                            model.addConstr(arrived,GRB.EQUAL,0,"No Airborne_Delay of "
                            +K[k]+" from "+J[jactu]+" to "+J[jnext]);
                        }
                        jactu = jnext;
                    }
                }
            }

            //Constraint "not more entry than capacity"
            System.out.println("Creating constraints no more entry than capacity");
            for (int l = 0 ;l<Sl.length;l++){
                for (int p= 0 ;p<np; p++){
                    GRBLinExpr entry_count_sum = GRBLinExpr();
                    for(int f=0; f<nf;f++){
                        for (int k:Kf.get(f)){
                            int j = Jktol[l][p][k];
                            int begin = max(Tfkj[k][j][0],P[p][0]);
                            int last = min(Tfkj[k][j][Tfkj[k][j].length -1],P[p][1]);
                            for (int t=begin ; t<last+1 ; t++){
                                entry_count_sum.addTerm(1,entry[k][j][t]);
                            } 
                        }
                    }
                    model.addConstr(entry_count_sum,GRB.LESS_EQUAL,cap[l][p],
                         " Limit capacity of sector l"+l+ " during period p"+p);
                }
            }

            long end_model_init_time = System.currentTimeMillis();
            long init_time = (end_model_init_time-start_model_init_time)/1000;
            System.out.println("The model was built in "+ 
            init_time + " s.");
        }
        catch(GRBException e){
            System.out.println("Error code: " + e.getErrorCode() + ". "+ 
            e.getMessage());
        }
        catch(IOException e){
            e.printStackTrace();
        }
    }


    public void dispose(){
        //dispose of model and environment
        model.dispose();
        env.dispose();
    }

    @Override
    public boolean solve(   String log_file,
                            int time_limit_seconds,
                            double heuristics_value,
                            int presolve_value,
                            int[][][] x_res,
                            int[] z_res,
                            int[][] u_res){
        boolean success = false;
        long begin_solve_time =  System.currentTimeMillis();

        // Solve model
        System.out.println("Solving");
        try{
            //TODO
            //Solver settings
            model.set(GRB.DoubleParam.Heuristics, heuristics_value);
            //model.set(GRB.IntParam.Threads,1);
            model.set(GRB.IntParam.Presolve, presolve_value);
            model.set(GRB.DoubleParam.TimeLimit, time_limit_seconds); //In seconds
            // TODO Callback?
            /*//Callback function
            System.err.println("Time-[ms] objbst objbnd");
            my_callback object_callback = new my_callback(model.getVars());
            model.setCallback(object_callback); */

            model.optimize();

            //object_callback.write_log_file(log_file); // Write the logs of the solver

            long end_solve_time =  System.currentTimeMillis();
            long solve_time = (end_solve_time-begin_solve_time)/1000;
            System.out.println("The model was solved in "+ 
                solve_time + " s.");

            if(model.get(GRB.IntAttr.Status) == GRB.Status.OPTIMAL) {
                success = true;
                System.out.println("Resolution successful");
                // Extract solution
                System.out.println("Extracting solution");
                x_res = new int[nk][nj][nt+1]; 
                for (int k = 0;k<nk;k++){
                    for (int j = 0; j <nj; j++){
                        for (int t = 0; t < nt+1; t++){
                            x_res[k][j][t] = x[k][j][t].get(GRB.DoubleAttr.X);
                        }
                    }
                }
                z_res = new int[nk];
                for (int k = 0;k<nk;k++){
                    z_res [k] = z[k].get(GRB.DoubleAttr.X);
                }
            }
            else {
                System.out.println("Resolution unsuccessful");
            }
        }
        catch (GRBException e) {
            System.out.println("Error code: " + e.getErrorCode() + ". " +
                    e.getMessage());
            e.printStackTrace();
        }
        return success;
    }


    @Override
    public void print_solution( double alph,
                                double gam,
                                double delt,
                                String[] F,
                                String[] K,  
                                String[] J,
                                double[][]r,    // r[k][j]
                                List<List<Integer>>  Kf, //Kf[f] is the list of indexes of flight f's trajectory in K
                                String[] C,    //array containing all the configurations
                                int [][] Ca, //index of configurations grouped by area control center 
                                int[][] Sl,  // Sl[l] gives an array listing the configurations containing the sector l
                                Integer[]  T, //T is the ordered set of time
                                int[][] P, // P = [...;[ind_begin;ind_end];...] periods of T
                                int[][][] Tfkj, // Tfkj[k][j] is an ordered array of time indexes where the flight following k can be in j the indexes are counted from 1
                                int[][]Jfki, // Jfki[k][i] is the ith elementary sector the flight following k passes through 
                                double[] D, // indexed by K
                                double[] E, // indexed by K
                                int [][] TC, // TC[k][j][jp] is the time needed (in number of time steps) for the flight following k to get from j to jp (sectors)
                                                // if k does not go from j to jp let TC[k][j][jp] be superior to T[].length
                                int [][] cap,  //cap[l][p] is the capacity of sector l during the period indexed p
                                int [][][] Jktol, // Jktol[l][p][k] represent the first elementary sector for k that function in l in time period p
                                int M,
                                int[][][] x_res,
                                int[] z_res,
                                int[][] u_res,
                                String output_filename) {

        int nt = T.length;
                                                
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(output_filename));

            // USED TRAJECTORIES
            System.out.println("Selected trajectories:");
            bw.write("Selected trajectories:\n");
            
            for (int f = 0; f<nf;f++){
                for (int k : Kf.get(f)){
                    if (z_res[k]==1){
                        int j = Jfki[k][0];
                        int t=1;
                        while (t<nt+1 && x_res[k][j][t]==x_res[k][j][t-1]){ // t<nt+1 is a safety but may not be needed
                            t++;
                        }
                        double delay =  (T[t-1]-r[k][j]);
                        double cost = gam*D[k]+E[k] + alph * delay;
                        System.out.println("Flight " + F[f] + " -> Trajectory " + K[k] + " -> Delay " + delay + " -> Cost " + cost);
                        bw.write("Flight " + F[f] + " -> Trajectory " + K[k] + " -> Delay " + delay + " -> Cost " + cost + "\n");
                    }
                }
            }

            bw.close();
        }
        catch (Exception e) {
            System.err.println("Problem in writing solution file! " + e);
        }
    }    
}