package data;
import java.io.*;
import java.util.*;
import java.util.stream.IntStream;
//import java.lang.Math.min;        //doesn't work for some reason

public class DataManager {
    public static void print_data_summary(String[] F,
                                          String[] K,
                                          String[] J,
                                          int[][] P,
                                          String[] C,
                                          Integer[] T, 
                                          int[][] Sl)
    {
        System.out.println("");
        System.out.println("Num. flights " + F.length);
        System.out.println("Num. trajectories " + K.length);
        System.out.println("Num. periods " + P.length);
        System.out.println("Num. times " + T.length);
        System.out.println("Num. sectors " + Sl.length);
        System.out.println("Num. elementary sectors " + J.length);
        //System.out.println("Num. clusters " + Cluster_set.size()); // what cluster?
        System.out.println("Num. configurations " + C.length);
        System.out.println("");
    }

    public static void initialise_data( double max_complx,
                                        int period_duration,
                                        String clusters_file,
                                        String mcx_file,
                                        String trajectory_cost_file,
                                        String configuration_file,
                                        String airspace_file,
                                        String sector_cost_file,
                                        double alph,
                                        double gam,
                                        double delt,
                                        String[] F,
                                        String[] K,  
                                        String[] J,
                                        double[][]r,    // r[k][j]
                                        List<List<Integer>> Kf, //Kf.get(f) is the list of indexes of flight f's trajectory in K
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
                                        int M,
                                        CapacitySelector selector){
        List<String> Se = new ArrayList<String>();
        List<String> S= new ArrayList<String>();
        Map<String,  List<String>> Eps = new HashMap<String,List<String>>();
        List<String> Cluster_set = new ArrayList<String>();
        Map<String,  Double> delta_cost = new HashMap<String,  Double>();
        Map<String, Integer> K_index;
        Map<String, Integer> J_index;
        String[] L;
        Map<String, Integer> L_index;

        long start_reading = System.currentTimeMillis();
        read_clusters_file(clusters_file,Cluster_set);
        read_trajectory_table_cost(trajectory_cost_file,delta_cost);
        read_configuration_file( configuration_file, Cluster_set, C, Ca, Sl,L);
        L_index = reverse(L);
        read_spc_file(airspace_file,S, Se, Eps,J);
        J_index =  reverse(J);

        Map<String,Double> theta_cluster = new HashMap<String,Double>();
        read_sector_cost_file(sector_cost_file,theta_cluster);
        read_data_from_mcx_file(mcx_file, max_complx,period_duration,Se,F,K,J,T,P,Kf,Tfkj,Jfki,TC,Jktol,K_index);
        cap = new int[L.length][P.length];
        Arrays.fill(cap,0);
    }

    public static Integer time_to_int(String time){
        String[] time_split = time.split(":");
        int res = Integer.parseInt(time_split[0])*60 + Integer.parseInt(time_split[1]) ;
        return res;
    }

    public static void read_ncap_file(String ncap_file,Integer[] T, int[][]P,Map<String, Integer> L_index, int[][]cap,CapacitySelector selector){
        try{

            System.out.println("Reading ncap file: " + ncap_file);
            BufferedReader reader = new BufferedReader(new FileReader(ncap_file));
            String line = "";
            String[] splitted_line;
            String sector;
            int l;
            int end;
            ArrayList<ArrayList<ArrayList<Integer>>> capacities = new ArrayList<ArrayList<ArrayList<Integer>>>(); 

            while ((line = reader.readLine()) != null) {
                splitted_line = line.split(";");
                sector = splitted_line[1];
                l = L_index.get(sector);
                end = time_to_int(splitted_line[3]);
                for (int p=0 ;p<P.length;p++){
                    int []plist = P[p];
                    if (T[p[0]]<= end){

                    }
                }

            }
            reader.close();
        }
        catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void read_spc_file(String spc_file,
                                     List<String> S,
                                     List<String> Se,
                                     Map<String,  List<String>> Eps,
                                     String[] J)
    {
        try {
            System.out.println("Reading spc file: " + spc_file);
            BufferedReader reader = new BufferedReader(new FileReader(spc_file));
            String line = "";
            String[] splitted_line;
            String airspace = "";
            String airspace_element ="";
            while ((line = reader.readLine()) != null) {
                splitted_line = line.split(";");
                if(Objects.equals(splitted_line[0], "A"))
                {
                    airspace = splitted_line[1];
                    if(S.contains(airspace))
                    {
                        Eps.put(airspace,new ArrayList<String>());
                    }
                }
                if(Objects.equals(splitted_line[0], "S"))
                {
                    airspace_element = splitted_line[1];
                    if(S.contains(airspace))
                    {
                        Eps.get(airspace).add(airspace_element);
                        //Se.add(airspace_element);
                        if(!Se.contains(airspace_element))
                        {
                            Se.add(airspace_element);
                        }
                    }
                }
            }
            reader.close();

            //Inlcude the elementary sectors
            for(String s:S)
            {
                if(!Eps.containsKey(s)) //In not in airspace file it means it is ES
                {
                    Eps.put(s,Arrays.asList(s));
                    if(!Se.contains(s))
                    {
                        Se.add(s);
                    }
                }
            }
            J = (String[]) Se.toArray();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static Map<String, Integer> reverse(String[] array){
        Map<String,Integer> res = new HashMap<String, Integer>();
        for (int i=0;i<array.length;i++){
            res.put(array[i], i);
        }
        return res;
    } 

    private static int[] transform(List<String> list, Map<String, Integer> indexes){
        int nl = list.size();
        int[] res = new int[nl];
        for (int ind =0;ind<nl;ind++){
            res[ind]= indexes.get(list.get(ind));
        }
        return res;
    }

    public static void read_configuration_file( String conf_file,
                                                List<String> Cluster_set,
                                                String[] C,    //array containing all the configurations
                                                int [][] Ca, //index of configurations grouped by area control center 
                                                int[][] Sl,  // Sl[l] gives an array listing the configurations containing the sector l                               
                                                String[] L)
    {
        try {
            System.out.println("Reading configuration file: " + conf_file);
            BufferedReader reader = new BufferedReader(new FileReader(conf_file));
            String line = "";
            String[] splitted_line;
            String ACC = "";
            String conf = "";
            String sector = "";

            Set<String> set_of_conf = new HashSet<String>();
            Set<String> set_of_sectors = new HashSet<String>();
            
            Map<String,  List<String>> C_s= new HashMap<String,  List<String>>();
            Map<String,  List<String>> conf_in_cluster= new HashMap<String,  List<String>>();
            Map<String,  String> cluster_conf= new HashMap<String,  String>();
            //Map<String,  List<String>> sectors_in_conf = new HashMap<String,  List<String>>();
            

            while ((line = reader.readLine()) != null) {
                splitted_line = line.split(";");
                ACC = splitted_line[0];
                conf = ACC + "_" + splitted_line[1]; //It ensures uniq ID
                sector = splitted_line[2];

                if(Cluster_set.contains(ACC))
                {
                    set_of_conf.add(conf);
                    set_of_sectors.add(sector); //It can be collapsed or elementary sector

                    if(!C_s.containsKey(sector))
                    {
                        C_s.put(sector,new ArrayList<String>());
                    }
                    C_s.get(sector).add(conf);

                    cluster_conf.put(conf,ACC);

                    /*if(!sectors_in_conf.containsKey(conf))
                    {
                        sectors_in_conf.put(conf,new ArrayList<String>());
                    }
                    sectors_in_conf.get(conf).add(sector);*/

                    if(!conf_in_cluster.containsKey(ACC))
                    {
                        conf_in_cluster.put(ACC,new ArrayList<String>());
                    }
                    if(!conf_in_cluster.get(ACC).contains(conf))
                    {
                        conf_in_cluster.get(ACC).add(conf);
                    }
                }

            }
            C = (String[]) set_of_conf.toArray();
            Map<String, Integer> Cind = reverse(C);
            int nacc=conf_in_cluster.size(); 
            Ca = new int[nacc][];
            int index=0;
            for (String key:conf_in_cluster.keySet()){
                Ca[index] = transform(conf_in_cluster.get(key), Cind);
                index++;
            }
            L = (String[]) set_of_sectors.toArray(); // if needed at some point, we just have to add it to the function's parameter
            Sl = new int[L.length][];
            for(int ind =0;ind<L.length;ind++){
                Sl[ind] = transform(C_s.get(L[ind]),Cind); 
            }
            
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
    public static void read_clusters_file(String clusters_file,List<String> Cluster_set)
    {
        try {
            System.out.println("Reading cluster file: " + clusters_file);
            BufferedReader reader = new BufferedReader(new FileReader(clusters_file));
            String line = "";
            String cluster = "";
            while ((line = reader.readLine()) != null) {
                cluster = line.trim();
                Cluster_set.add(cluster);
            }
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static void read_data_from_mcx_file(String mcx_file,
                                                double max_complx,
                                                int period_duration,
                                                List<String> Se,
                                                String[] F,
                                                String[] K,
                                                String[] J,
                                                Integer[] T,
                                                int[][] P,
                                                List<List<Integer>> Kf,
                                                int[][][] Tfkj,
                                                int[][]Jfki,
                                                int [][] TC,
                                                int [][][] Jktol,
                                                Map<String, Integer> K_index)
    {
        try {
            System.out.println("Reading mcx file: " + mcx_file);
            BufferedReader reader = new BufferedReader(new FileReader(mcx_file));
            String line = "";
            String[] splitted_line;
            String trajectory_i, flight_i, time, elem_sector;
            String trajectory_j, flight_j;
            double complx;
            Set<String> F_set = new HashSet<String>();
            Set<String> K_set = new HashSet<String>();
            Set<String> T_set = new HashSet<String>();
            Map<String,  List<String>> Kf_map = new HashMap<String,  List<String>>();

            while ((line = reader.readLine()) != null) {
                splitted_line = line.split(",");
                trajectory_i = splitted_line[0];
                flight_i = trajectory_i.split("_")[0];
                time = splitted_line[1];
                elem_sector = splitted_line[2];
                if(Se.contains(elem_sector)) // Only if it is the cluster zone
                {
                    F_set.add(flight_i);
                    K_set.add(trajectory_i);
                    T_set.add(time);

                    if(!Kf_map.containsKey(flight_i))
                    {
                        Kf_map.put(flight_i,new ArrayList<String>());
                    }
                    if(!Kf_map.get(flight_i).contains(trajectory_i))
                    {
                        Kf_map.get(flight_i).add(trajectory_i);
                    }

                    if(splitted_line.length>3)
                    {
                        for(int i = 3; i<splitted_line.length-1;i=i+2)
                        {
                            trajectory_j = splitted_line[i];
                            //flight_j = trajectory_j.split("_")[0];

                            // Not considered. One trajectory outise the cluste may influence inside, but no regulation applies there
//                            K_set.add(trajectory_j); // Is it needed? If appear here it should appear as trajectory_i
//                            F_set.add(flight_j);     // Is it needed? If appear here it should appear as trajectory_i
//                            if(!Kf.containsKey(flight_j))
//                            {
//                                Kf.put(flight_j,new ArrayList<String>());
//                            }
//                            if(!Kf.get(flight_j).contains(trajectory_j))
//                            {
//                                Kf.get(flight_j).add(trajectory_j);
//                            }

                            complx = Double.parseDouble(splitted_line[i+1]);
                            if(complx > max_complx)
                            {
                                complx = max_complx;
                            }
                        }
                    }
                }
            }

            List<Integer> t_int_vector = new ArrayList<Integer>();
            for(String t:T_set)
            {
                t_int_vector.add(Integer.parseInt(t));
            }
            Collections.sort(t_int_vector);
            //int period_duration = 10*60; // 10 minutes
            int num_times_in_period = period_duration*60/20;
            //num_times_in_period = 1; //TOY EXAMPLE
            int min_period = t_int_vector.get(0)/num_times_in_period+1;
            int max_period = t_int_vector.get(t_int_vector.size()-1)/num_times_in_period+1;
            System.out.println("Min. period: p"+min_period);
            System.out.println("Max. period: p"+max_period);
            P = new int[max_period-min_period+1][2];
            T = (Integer[])t_int_vector.toArray();
            int maxi = t_int_vector.size();
            for(int p_int = 0; p_int <= max_period-min_period; p_int++)
            {
                P[p_int][0] = p_int*num_times_in_period;
                P[p_int][1] = min((p_int+1)*num_times_in_period,maxi)-1;
            }

            F = (String[]) F_set.toArray();
            K = (String[]) K_set.toArray();

            K_index =  reverse(K);
            Kf = new ArrayList<List<Integer>>();
            for (int f=0;f<F.length;f++){
                List<Integer> Kfks = new ArrayList<Integer>();
                List<String> interm = Kf_map.get(F[f]);
                for(String k: interm){
                    Kfks.add(K_index.get(k));
                }
                Kf.add(f, Kfks);
            }


            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static int min(int i, int j) {
        if (i<j){return i;}
        else{return j;}
    }

}
