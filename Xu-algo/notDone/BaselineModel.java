import gurobi.*;

import java.io.*;
import java.util.*;

public class BaselineModel implements Model{

    public static String model_name () {return "Baseline Model"}

    protected static boolean solveModel(List<String> F,     // list of flight
                                        List<String> K,     //
                                        List<String> Cluster_set,
                                        List<String> C,
                                        List<String> S,
                                        List<String> Se,    // I suppose this is the list of elementary sectors
                                        List<String> T,     // I suppose this is the list of time steps
                                        List<String> P,
                                        Map<String,  List<String>> times_in_p,
                                        Map<String,  List<Integer>> ind_times_in_p,
                                        Map<String,  List<String>> Kf,
                                        Map<String,  List<String>> C_s,
                                        Map<String,  List<String>> conf_in_cluster,
                                        Map<String,  List<String>> sectors_in_conf,
                                        Map<String,  List<String>> Eps,
                                        Map<String,  Double> delta_cost,
                                        Map<String,  Map<String,Double>> theta,
                                        double phi,
                                        Map<String, Map<String,Double>> threshold,
                                        int M,
                                        int M_h,
                                        Map<String,Map<String,Map<String,Double>>> C_kk,
                                        Map<String,Map<String,Map<String,Integer>>> B,
                                        Map<String,  Double> z_var,
                                        Map<String,  Map<String,Double>> x_var,
                                        Map<String,  Map<String,Double>> y_var,
                                        Map<String,Map<String,Map<String,Double>>> Ct_var,
                                        Map<String,  Map<String,Double>> q_var,
                                        String output_filename){
        boolean success = false;
        int nf = F.length();
        int nt = T.size();
        int nj = Se.size();
        try{
            //Model
            long start_model_time = System.currentTimeMillis();
            GRBEnv env = new GRBEnv();
            GRBModel model = new GRBModel(env);
            model.set(GRB.StringAttr.ModelName, "baseline model");

            // Create decision variable for delay selection
            System.out.println("Creating decision variable x");
            GRBVar[][][] x = new GRBVar[nf][nt+1][nj];                       // +1 because we consider the time step 0
            for (int f = 0; f < nf; f++ ){
                for (int t = 0; t < nt+1; t++){
                    for (int j = 0; j <nj; j++){
                        x[f][t][j] = model.addVar(0,1,0,GRB.BINARY,
                            "x_f"+((String)f)+"_t"+((String)t)+"_j"+((String)j));
                    }
                }
            }

            //Objective function
            //TODO
            System.out.println("Creating objective function");
            GRBLinExpr obj = new GRBLinExpr();
                //cost of delay
            

            //Constraint "cannot arrive before possible"
            //TODO
            //Constraint "must arrive before last possible time"
            //TODO
            //Constraint "if we previously arrived, it id true in the future"
            //TODO
            //Constraint "no airborne delay"
            //TODO
            //Constraint "not more entry than capacity"
            //TODO

            //Solve model
            //TODO

            //dispose of model and environment
            long end_model_time = System.currentTimeMillis();
            model.dispose();
            env.dispose();
            System.out.println("The model was running for "+ 
                (String)(end_model_time-start_model_time) + " ms.");
        }
        catch(GRBException e){
            System.out.println("Erro code: " + e.getErrorCode() + ". "+ 
                e.getMessage());
        }
        catch(IOException e){
            e.printStackTrace();
        }

    }







}