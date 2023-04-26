package notDone;
import java.io.*;
import java.util.*;
import java.util.stream.IntStream;

public class DataManagerOld {

    public static void print_data_summary(String[] F,
                                          String[] K,
                                          String[] J,
                                          int[][] P,
                                          String[] C,
                                          double[] T, 
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
/*
    public static void initialise_data(double max_complx,
                                       int period_duration,
                                       String clusters_file,
                                       String mcx_file,
                                       String trajectory_cost_file,
                                       String configuration_file,
                                       String airspace_file,
                                       String sector_cost_file,
                                       List<String> F,
                                       List<String> K,
                                       List<String> Cluster_set,
                                       List<String> C,
                                       List<String> S,
                                       List<String> Se,
                                       List<String> T,
                                       List<String> P,
                                       Map<String,  List<String>> times_in_p,
                                       Map<String,  List<Integer>> ind_times_in_p,
                                       Map<String,  List<String>> Kf,
                                       Map<String,  List<String>> Eps,
                                       Map<String,  List<String>> C_s,
                                       Map<String,  List<String>> conf_in_cluster,
                                       Map<String,  String> cluster_conf,
                                       Map<String,  List<String>> sectors_in_conf,
                                       Map<String,  Double> delta_cost,
                                       double threshold_value,
                                       Map<String, Map<String,Double>> threshold,
                                       Map<String,Map<String,Map<String,Double>>> C_kk,
                                       Map<String,Map<String,Map<String,Integer>>> B,
                                       Map<String,  Map<String,Double>> theta_period)
    {
    //        Cluster_set.add("EDUUUTAC");
    //        Cluster_set.add("EDUUUTAE");

        long start_reading = System.currentTimeMillis();
        read_clusters_file(clusters_file,Cluster_set);
        read_trajectory_table_cost(trajectory_cost_file,delta_cost);
        read_configuration_file(configuration_file, C, S, Cluster_set, C_s, conf_in_cluster, cluster_conf, sectors_in_conf);

        read_spc_file(airspace_file,S, Se, Eps);

        Map<String,Double> theta_cluster = new HashMap<String,Double>();
        read_sector_cost_file(sector_cost_file,theta_cluster);
        read_data_from_mcx_file(mcx_file, max_complx,period_duration,F,K,Se,T,P,times_in_p,ind_times_in_p,Kf,C_kk,B);

        for(String t:T)
        {
            threshold.put(t,new HashMap<String,Double>());
            for(String s:S)
            {
                //threshold.get(t).put(s,1300.0);
                //threshold.get(t).put(s,2.0); //TOY EXAMPLE
                threshold.get(t).put(s,threshold_value);
            }
        }

        // Cost of opening sector s during the time t - Index: {t,s}
        // Cost of opening sector s during the time t - Index: {t,s}
        Map<String,  Map<String,Double>> theta = new HashMap<String,  Map<String,Double>>();
        //theta_value = 100.0;
        for(String t:T)
        {
            theta.put(t,new HashMap<String,Double>());
            for(String s:S)
            {
                //theta.get(t).put(s,theta_value);


                String conf = C_s.get(s).get(0); //I consider only one conf
                String clust = cluster_conf.get(conf);
                // 1t last 20 seconds -> division by 180 the hourly cost
                theta.get(t).put(s,theta_cluster.get(clust)/180);
            }
        }
        for(String p:P)
        {
            theta_period.put(p,new HashMap<String,Double>());
            for(String s:S)
            {
                theta_period.get(p).put(s,0.0);
                for(String t:times_in_p.get(p))
                {
                    theta_period.get(p).put(s,theta_period.get(p).get(s)+theta.get(t).get(s));
                }
            }
        }

        long end_reading = System.currentTimeMillis();
        long elapsedTime_reading = (end_reading - start_reading)/1000; //Seconds
        System.out.println();
        System.out.println("Time for initialise data: " + elapsedTime_reading + " s");

    }

    private static void read_sector_cost_file(String sector_cost_file,Map<String,Double> theta_cluster)
    {
        try {
            System.out.println("Reading sector cost file: " + sector_cost_file);
            BufferedReader reader = new BufferedReader(new FileReader(sector_cost_file));
            String[] splitted_line;
            String line = "";
            String cluster = "";
            double sector_cost;
            reader.readLine(); //Header
            while ((line = reader.readLine()) != null) {
                splitted_line = line.split(";");
                cluster = splitted_line[0];
                sector_cost = Double.parseDouble(splitted_line[6]);
                theta_cluster.put(cluster,sector_cost);
            }

            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    

    private static void read_data_from_mcx_file(String mcx_file,double max_complx,int period_duration,
                                               List<String> F,
                                               List<String> K,
                                               List<String> Se,
                                               List<String> T,
                                               List<String> P,
                                               Map<String,  List<String>> times_in_p,
                                               Map<String,  List<Integer>> ind_times_in_p,
                                               Map<String,  List<String>> Kf,
                                               Map<String,Map<String,Map<String,Double>>> C_kk,
                                               Map<String,Map<String,Map<String,Integer>>> B)
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

                    if(!Kf.containsKey(flight_i))
                    {
                        Kf.put(flight_i,new ArrayList<String>());
                    }
                    if(!Kf.get(flight_i).contains(trajectory_i))
                    {
                        Kf.get(flight_i).add(trajectory_i);
                    }

                    if(!B.containsKey(trajectory_i))
                    {
                        B.put(trajectory_i,new HashMap<String,Map<String,Integer>>());
                    }
                    if(!B.get(trajectory_i).containsKey(time))
                    {
                        B.get(trajectory_i).put(time, new HashMap<String,Integer>());
                    }
                    B.get(trajectory_i).get(time).put(elem_sector,1);

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

                            if(!C_kk.containsKey(trajectory_i))
                            {
                                C_kk.put(trajectory_i,new HashMap<String,Map<String,Double>>());
                            }
                            if(!C_kk.get(trajectory_i).containsKey(trajectory_j))
                            {
                                C_kk.get(trajectory_i).put(trajectory_j,new HashMap<String,Double>());
                            }
                            C_kk.get(trajectory_i).get(trajectory_j).put(time,complx);
                        }
                    }
                }


            }

            List<Integer> t_int_vector = new ArrayList<>();
            for(String t:T_set)
            {
                t_int_vector.add(Integer.parseInt(t));
            }
            //int period_duration = 10*60; // 10 minutes
            int num_times_in_period = period_duration*60/20;
            //num_times_in_period = 1; //TOY EXAMPLE
            int min_period = Collections.min(t_int_vector)/num_times_in_period+1;
            int max_period = Collections.max(t_int_vector)/num_times_in_period+1;
            System.out.println("Min. period: p"+min_period);
            System.out.println("Max. period: p"+max_period);
            for(int p_int = min_period; p_int <= max_period; p_int++)
            {
                String p = "p"+p_int;
                P.add(p);
                times_in_p.put(p, new ArrayList<String>());
                for(int t_int = (p_int-1)*num_times_in_period; t_int < p_int*num_times_in_period; t_int++)
                {
                    String t = String.valueOf(t_int);
                    times_in_p.get(p).add(t);
                    T.add(t);
                }
            }

            // List of index (in T) of the times in period p. Key: {p}
            for(String p:times_in_p.keySet())
            {
                ind_times_in_p.put(p,new ArrayList<Integer>());
                for(String t:times_in_p.get(p))
                {
                    int t_ind = T.indexOf(t);
                    ind_times_in_p.get(p).add(t_ind);
                }
            }

            F.addAll(F_set);
            K.addAll(K_set);

            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static void read_spc_file(String spc_file,
                                     List<String> S,
                                     List<String> Se,
                                     Map<String,  List<String>> Eps)
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
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static void read_clusters_file(String clusters_file,List<String> Cluster_set)
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

    private static void read_configuration_file(String conf_file,
                                               List<String> C,
                                               List<String> S,
                                               List<String> Cluster_set,
                                               Map<String,  List<String>> C_s,
                                               Map<String,  List<String>> conf_in_cluster,
                                               Map<String,  String> cluster_conf,
                                               Map<String,  List<String>> sectors_in_conf)
    {
        try {
            System.out.println("Reading configuration file: " + conf_file);
            BufferedReader reader = new BufferedReader(new FileReader(conf_file));
            String line = "";
            String[] splitted_line;
            String ACC = "";
            String conf = "";
            String sector = "";

            Set<String> set_of_configurations = new HashSet<String>();
            Set<String> set_of_sectors = new HashSet<String>();

            while ((line = reader.readLine()) != null) {
                splitted_line = line.split(";");
                ACC = splitted_line[0];
                conf = ACC + "_" + splitted_line[1]; //It ensures uniq ID
                sector = splitted_line[2];

                if(Cluster_set.contains(ACC))
                {
                    set_of_configurations.add(conf);
                    set_of_sectors.add(sector); //It can be collapsed or elementary sector

                    if(!C_s.containsKey(sector))
                    {
                        C_s.put(sector,new ArrayList<String>());
                    }
                    C_s.get(sector).add(conf);

                    cluster_conf.put(conf,ACC);

                    if(!sectors_in_conf.containsKey(conf))
                    {
                        sectors_in_conf.put(conf,new ArrayList<String>());
                    }
                    sectors_in_conf.get(conf).add(sector);

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
            C.addAll(set_of_configurations);
            S.addAll(set_of_sectors);
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static void read_trajectory_table_cost(String filename,Map<String,  Double> delta_cost)
    {
        try {
            System.out.println("Reading trajectory cost table: " + filename);
            BufferedReader reader = new BufferedReader(new FileReader(filename));
            String line = "";
            String[] splitted_line;
            boolean first_line = true;

            while ((line = reader.readLine()) != null) {
                if(first_line)
                {
                    first_line = false;
                }
                else
                {
                    splitted_line = line.split(";");
                    String trajectory_id = splitted_line[0];
                    double delta_total_cost = Double.parseDouble(splitted_line[19]);
                    delta_cost.put(trajectory_id,delta_total_cost);
                }
            }
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
*/
}
