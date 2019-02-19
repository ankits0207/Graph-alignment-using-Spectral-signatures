import java.io.IOException;
import java.util.*;

import org.jblas.ComplexDouble;
import org.jblas.ComplexDoubleMatrix;
import org.jblas.DoubleMatrix;
import org.jblas.Eigen;

public class SeedPhase {
	public ArrayList<Score> get_dist_among_all_pairs(Graph g1, Graph g2, int k, double alpha, String eval_file_name, double default_eval) throws IOException
	{
		Seed seed = new Seed();
		seed.score= 100000;
		
		ArrayList<Score> al_of_score = new ArrayList<Score>();
		
		BLAST_evalues all_align = new BLAST_evalues();
		ArrayList<Hit> al_of_hit = all_align.read_Evalues(eval_file_name);
		
		int i = 0;
		int total = g1.al_of_nodes.size()*g2.al_of_nodes.size();
		double topo_dist, seq_sim;
		
		for(MyNode mn1 : g1.al_of_nodes)
		{
			for(MyNode mn2 : g2.al_of_nodes)
			{
				System.out.println("********************" + i + "********************" + total);
				Score s = new Score();
				s.n1 = mn1;
				s.n2 = mn2;
				if(alpha == 1.0)
				{
					s.score = get_k_hop_topo_dist(g1, g2, mn1, mn2, k);
				}
				else
				{
					topo_dist = get_k_hop_topo_dist(g1, g2, mn1, mn2, k);
					seq_sim = get_sim_score(mn1, mn2, al_of_hit, default_eval);
					s.score = get_d_alpha(alpha, topo_dist, seq_sim);
				}
				al_of_score.add(s);
				i++;
			}
		}
		
		return al_of_score;
	}
	
	public double get_sim_score(MyNode mn1, MyNode mn2, ArrayList<Hit> al_of_hit, double default_eval)
	{
		double sim_score = default_eval;
		for(Hit h : al_of_hit)
		{
			if(h.Node1.equals(mn1.Id) && h.Node2.equals(mn2.Id))
			{
				sim_score  = h.evalue;
			}
		}
		return sim_score;
	}
	
	public double get_topological_distance(Graph g1, Graph g2, HashMap<String, Integer> hm1, HashMap<String, Integer> hm2){
		Get_Matrices gm = new Get_Matrices();
		Utility u = new Utility();
		
		int[][] adj_matrix_1 = gm.get_adjacency_matrix(g1, hm1);
		int[][] adj_matrix_2 = gm.get_adjacency_matrix(g2, hm2);
		
		double[][] normalized_laplacian_1 = gm.get_normalized_laplacian_matrix(adj_matrix_1, u.get_degree_map(adj_matrix_1, hm1));
		double[][] normalized_laplacian_2 = gm.get_normalized_laplacian_matrix(adj_matrix_2, u.get_degree_map(adj_matrix_2, hm2));
		
		DoubleMatrix dm1 = new DoubleMatrix(normalized_laplacian_1);
		DoubleMatrix dm2 = new DoubleMatrix(normalized_laplacian_2);
		
		ComplexDoubleMatrix eigen_values_1 = Eigen.eigenvalues(dm1);
		ComplexDoubleMatrix eigen_values_2 = Eigen.eigenvalues(dm2);
		
		double[] eigen_val_decomp_1 = new double[eigen_values_1.toArray().length];
		double[] eigen_val_decomp_2 = new double[eigen_values_2.toArray().length];
		
		int idx = 0;
		for(ComplexDouble cde1 : eigen_values_1.toArray())
		{
			eigen_val_decomp_1[idx] = cde1.abs();
			idx++;
		}
		
		idx = 0;
		for(ComplexDouble cde2 : eigen_values_2.toArray())
		{
			eigen_val_decomp_2[idx] = cde2.abs();
			idx++;
		}
		
		//System.out.println("Binning Eigen Values...");
		int n = 100;
		double bin_size = 2.0/n;
		int freq;
		HashMap<Integer, Integer> hm_from_ev_1 = new HashMap<Integer, Integer>();
		for(int i = 1; i < n; i++)
		{
			hm_from_ev_1.put(i, 0);
		}
		HashMap<Integer, Integer> hm_from_ev_2 = new HashMap<Integer, Integer>();
		for(int i = 1; i < n; i++)
		{
			hm_from_ev_2.put(i, 0);
		}
		
		for(Double d : eigen_val_decomp_1)
		{
			for(Integer i : hm_from_ev_1.keySet())
			{
				if((i-1)*bin_size <= d && d < i*bin_size)
				{
					freq = hm_from_ev_1.get(i);
					freq++;
					hm_from_ev_1.put(i, freq);
				}
				else if(d == i*bin_size)
				{
					freq = hm_from_ev_1.get(i);
					freq++;
					hm_from_ev_1.put(i, freq);					
				}
			}
		}
		for(Double d : eigen_val_decomp_2)
		{
			for(Integer i : hm_from_ev_2.keySet())
			{
				if((i-1)*bin_size <= d && d < i*bin_size)
				{
					freq = hm_from_ev_2.get(i);
					freq++;
					hm_from_ev_2.put(i, freq);
				}
				else if(d == i*bin_size)
				{
					freq = hm_from_ev_2.get(i);
					freq++;
					hm_from_ev_2.put(i, freq);					
				}
			}
		}
		
		//System.out.println("Transferring bin data to array...");
		double[] da1 = new double[n];
		double[] da2 = new double[n];
		int j = 0;
		for(Integer i : hm_from_ev_1.keySet())
		{
			da1[j] = hm_from_ev_1.get(i);
			j++;
		}
		j=0;
		for(Integer i : hm_from_ev_2.keySet())
		{
			da2[j] = hm_from_ev_2.get(i);
			j++;
		}
		return jensenShannonDivergence(da1, da2);
	}

	public static double jensenShannonDivergence(double[] p1, double[] p2) {
		//System.out.println("Finding Jensen Shannon divergence");
		assert(p1.length == p2.length);
	    double[] average = new double[p1.length];
	    for (int i = 0; i < p1.length; ++i) {
	    	average[i] += (p1[i] + p2[i])/2;
	    }
	    return (klDivergence(p1, average) + klDivergence(p2, average))/2;
	}

	public static final double log2 = Math.log(2);
	   
	public static double klDivergence(double[] p1, double[] p2) {
		double klDiv = 0.0;
		for (int i = 0; i < p1.length; ++i) {
			if (p1[i] == 0) { continue; }
			if (p2[i] == 0.0) { continue; }
			klDiv += p1[i] * Math.log( p1[i] / p2[i] );
		}
    return klDiv / log2;
    }
	
	public double get_k_hop_topo_dist(Graph g1, Graph g2, MyNode n1, MyNode n2, int k)
	{
		double d_topo = 0.0;
		Utility u = new Utility();		
		
		for(int i = 1; i <= k; i++)
		{
			Graph g_one = u.get_k_hop_neighbourhood(g1, n1, i);
			Graph g_two = u.get_k_hop_neighbourhood(g2, n2, i);
			HashMap<String, Integer> hash_map_1;
			HashMap<String, Integer> hash_map_2;
			hash_map_1 = u.get_HashMap(g_one);
			hash_map_2 = u.get_HashMap(g_two);
			d_topo += get_topological_distance(g_one, g_two, hash_map_1, hash_map_2);
		}
		return d_topo;
	}
	
	public double get_d_alpha(double alpha, double topo, double sim) throws IOException
	{
		double d_alpha;
		d_alpha = alpha*(topo) + (1.0 - alpha)*sim;
		return d_alpha;
	}
}
class Seed{
	MyNode n1;
	MyNode n2;
	double score;
}

class Score{
	MyNode n1;
	MyNode n2;
	double score;
	int flag;
}
