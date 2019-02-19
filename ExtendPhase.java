import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Set;

import org.jblas.*;

public class ExtendPhase {
	public ArrayList<Alignment> extend_alignment(Seed s, Graph g1, Graph g2, ArrayList<Score> al_of_score, int k, int score_idx){
		ArrayList<Alignment> al_of_a = new ArrayList<Alignment>();
		Queue<Seed> q = new LinkedList<Seed>();
		q.add(s);
		
		ArrayList<MyNode> al_of_visited_nodes_1 = new ArrayList<MyNode>();
		ArrayList<MyNode> al_of_visited_nodes_2 = new ArrayList<MyNode>();
		
		int number_of_nodes_g1 = g1.al_of_nodes.size();
		int number_of_nodes_g2 = g2.al_of_nodes.size();
		
		while(al_of_visited_nodes_1.size() < number_of_nodes_g1 && al_of_visited_nodes_2.size() < number_of_nodes_g2)
		{
			if(q.isEmpty())
			{
				score_idx++;
				Score score = al_of_score.get(score_idx);
				
				while(check_presence(al_of_visited_nodes_1, score.n1) != 0 || check_presence(al_of_visited_nodes_2, score.n2) != 0)
				{
					score_idx++;
					score = al_of_score.get(score_idx);	
				}
				
				Seed seed = new Seed();
				seed.n1 = score.n1;
				seed.n2 = score.n2;
				seed.score = score.score;
				q.add(seed);
			}
			System.out.println("*****SIZE OF ARRAYLIST OF ALIGNMENT*****" + al_of_a.size());
			System.out.println("*****FIRST*****" + al_of_visited_nodes_1.size());
			System.out.println("*****SECOND*****" + al_of_visited_nodes_2.size());
			
			Seed head_seed = q.remove();
			
			Utility u = new Utility();
			Graph subg1 = u.get_k_hop_neighbourhood(g1, head_seed.n1, 1);
			Graph subg2 = u.get_k_hop_neighbourhood(g2, head_seed.n2, 1);
			
			for(MyNode node : al_of_visited_nodes_1)
			{
				if(check_presence(subg1.al_of_nodes, node) == 1)
				{
					subg1.al_of_nodes.remove(node);
				}
			}
			
			for(MyNode node : al_of_visited_nodes_2)
			{
				if(check_presence(subg2.al_of_nodes, node) == 1)
				{
					subg2.al_of_nodes.remove(node);
				}
			}
			
			int size_of_q_mat = subg1.al_of_nodes.size()*subg2.al_of_nodes.size();
			
			if(size_of_q_mat > 0)
			{
				ArrayList<Pair> al_of_p = new ArrayList<Pair>();
				for(MyNode mn1 : subg1.al_of_nodes)
				{
					for(MyNode mn2 : subg2.al_of_nodes)
					{
						Pair p = new Pair();
						p.elt1 = mn1;
						p.elt2 = mn2;
						al_of_p.add(p);
					}
				}
				
				//Creating q matrix
				double[][] q_mat = new double[size_of_q_mat][size_of_q_mat];
				
				for(int i = 0; i < size_of_q_mat; i++)
				{
					for(int j = 0; j < size_of_q_mat; j++)
					{
						if(i == j)
						{
							Pair p = al_of_p.get(i);
							q_mat[i][j] = 1.0 - get_score(p.elt1, p.elt2, al_of_score);
						}
						else
						{
							Pair p1 = al_of_p.get(i);
							Pair p2 = al_of_p.get(j);
							q_mat[i][j] = get_c_val(p1, p2, al_of_score, g1, g2, k);
						}
					}
				}
				
				DoubleMatrix dm = new DoubleMatrix(q_mat);
				ComplexDoubleMatrix eigen_values = org.jblas.Eigen.eigenvalues(dm);

				
				double max_eigen_val = -100;
				int i = 0;
				int idx = 0;
				for(ComplexDouble eigen_val : eigen_values.toArray())
				{
					if(eigen_val.abs() > max_eigen_val)
					{
						max_eigen_val = eigen_val.abs();
						idx = i;
					}
					i++;
				}
				
				List<Double> principal_eigen_vector = getPrincipalEigenVector(dm, idx);
				int[][] solution_vector = new int[size_of_q_mat][1];
				
				ArrayList<Pair> temp_list = new ArrayList<Pair>();
				for(Pair p : al_of_p)
				{
					temp_list.add(p);
				}
				
				//Getting solution vector
				solution_vector = get_solution_vector(principal_eigen_vector, solution_vector, temp_list);
				
				//Assigning Nodes
				for(int j = 0; j < size_of_q_mat; j++)
				{
					Pair p = al_of_p.get(j);
					if(solution_vector[j][0] == 1)
					{
						Seed probable_seed = new Seed();
						probable_seed.n1 = p.elt1;
						probable_seed.n2 = p.elt2;
						q.add(probable_seed);
						
						Alignment a = new Alignment();
						a.node1 = p.elt1.Id;
						a.node2 = p.elt2.Id;

						int flg = 0;
						for(Alignment align : al_of_a)
						{
							if(align.node1.equals(a.node1) || align.node2.equals(a.node2))
							{
								flg = 1;
								break;									
							}
						}
						
						if(flg == 0)
						{
							al_of_a.add(a);
						}
					}
					al_of_visited_nodes_1.add(p.elt1);
					al_of_visited_nodes_2.add(p.elt2);
				}				
			}
			Set<MyNode> myset1 = new HashSet<MyNode>();
			myset1.addAll(al_of_visited_nodes_1);
			al_of_visited_nodes_1.clear();
			al_of_visited_nodes_1.addAll(myset1);
			
			Set<MyNode> myset2 = new HashSet<MyNode>();
			myset2.addAll(al_of_visited_nodes_2);
			al_of_visited_nodes_2.clear();
			al_of_visited_nodes_2.addAll(myset2);
		}
		return al_of_a;
	}
	
	public double get_score(MyNode mn1, MyNode mn2, ArrayList<Score> al_of_s)
	{
		double score = 0.0;
		for(Score s : al_of_s)
		{
			if(s.n1.Id.equals(mn1.Id) && s.n2.Id.equals(mn2.Id))
			{
				score = s.score;
				break;
			}
		}
		return score;		
	}

	public double get_c_val(Pair p1, Pair p2, ArrayList<Score> al_of_s, Graph g1, Graph g2, int k)
	{
		double c_val = 0.0;
		SeedPhase s1 = new SeedPhase();
		double d1 = s1.get_k_hop_topo_dist(g1, g2, p1.elt1, p2.elt1, k);
		double d2 = s1.get_k_hop_topo_dist(g1, g2, p1.elt2, p2.elt2, k);
		
		//double d1 = get_score(p1.elt1, p2.elt1, al_of_s);
		//double d2 = get_score(p1.elt2, p2.elt2, al_of_s);
		
		double num = -Math.abs(d1 - d2);
		double den = d1 + d2;
		if(den == 0)
		{
			return 0.0;		
		}
		else
		{
			c_val = Math.exp(num/den);
			return c_val;
		}	
	}
	
	public List<Double> getPrincipalEigenVector(DoubleMatrix dm, int idx)
	{
		ComplexDoubleMatrix eigenVectors = org.jblas.Eigen.eigenvectors(dm)[0];
		return getEigenVector(eigenVectors, idx);
	}
	
	public List<Double> getEigenVector(ComplexDoubleMatrix eigenvector, int columnId) {
	    ComplexDoubleMatrix column = eigenvector.getColumn(columnId);
	 
	    List<Double> values = new ArrayList<Double>();
	    for (ComplexDouble value : column.toArray()) {
	        values.add(value.abs());
	    }
	    return values;
	}
	
	public int[][] get_solution_vector(List<Double> eigen_vector_list, int[][] sol_vec, ArrayList<Pair> L){
		
		ArrayList<CustomPair3> al_of_cp3 = new ArrayList<CustomPair3>();
		for(int i = 0; i < eigen_vector_list.size(); i++)
		{
			CustomPair3 cp3 = new CustomPair3();
			cp3.eval = eigen_vector_list.get(i);
			cp3.p = L.get(i);
			al_of_cp3.add(cp3);			
		}
				
		int index = 0;
		int controller = 0;
		
		while(get_max(al_of_cp3).max != 0.0)
		{
			index = get_max(al_of_cp3).idx;
			
			sol_vec[index][0] = 1;
			controller = 0;
			Pair p = al_of_cp3.get(index).p;
			
			al_of_cp3.get(index).flag = 1;
			
			for(int j = 0; j < al_of_cp3.size(); j++)
			{
				Pair pair = al_of_cp3.get(j).p;
				if(al_of_cp3.get(j).flag == 0 && (pair.elt1.Id.equals(p.elt1.Id) || pair.elt2.Id.equals(p.elt2.Id)))
				{
					al_of_cp3.get(j).flag = 1;
				}
			}

			for(CustomPair3 cp3 : al_of_cp3)
			{
				if(cp3.flag == 0)
				{
					controller = 1;
				}
			}
			if(controller == 0)
			{
				break;
			}
		}
		return sol_vec;
	}
	
	public CustomPair1 get_max(ArrayList<CustomPair3> cp3_list){
		int idx = 0;
		int i = 0;
		Double max = -1000000.0;
		for(CustomPair3 cp3obj : cp3_list)
		{
			if(cp3obj.flag == 0 && cp3obj.eval > max)
			{
				max = cp3obj.eval;
				idx = i;
			}
			i++;
		}
		CustomPair1 cp = new CustomPair1();
		cp.max = max;
		cp.idx = idx;
		return cp;
	}
	
	public double get_ev(List<Double> al_of_d, int idx)
	{
		return al_of_d.get(idx);
	}
	
	public int check_presence(ArrayList<MyNode> al_of_n, MyNode mn)
	{
		for(MyNode n : al_of_n)
		{
			if(n.Id.equals(mn.Id))
			{
				return 1;
			}
		}
		return 0;
	}
}

class Alignment{
	String node1;
	String node2;
}

class Pair{
	MyNode elt1;
	MyNode elt2;
}

class CustomPair1
{
	int idx;
	double max;
}

class CustomPair2
{
	Pair p;
	double ev;
}

class CustomPair3
{
Pair p;
double eval;
int flag;
}
