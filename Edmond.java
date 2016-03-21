import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Scanner;

/*
 * Implementation of Edmond's Branching algorithm
 * Group Information: G13
 * Name : Diksha Sharma
 * Name : Sharol Clerit Pereira
*/

//Reference: https://www.cs.princeton.edu/courses/archive/spring13/cos423/lectures/04DemoEdmondsBranching.pdf
public class Edmond 
{
	public static Graph g; //Input graph
	public static int iOriginalVerticesCount = 0; //Original number of vertices in the graph
	public static long MSTWeight = 0; //MST weight of the graph
	public static ArrayList<ArrayList<Edge>> contractedCycles = new ArrayList<ArrayList<Edge>>();//List of List of edges that were contracted each time
	public static ArrayList<ArrayList<Vertex>> contractedVertices = new ArrayList<ArrayList<Vertex>>();//List of List of vertices that were contracted each time
	public static ArrayList<Edge> MST = new ArrayList<Edge>(); //Edges that form MST
	static int iContractedCycles = 0;//Number of times we contracted the cycles for the input graph
	private static int phase = 0;
	private static long startTime, endTime, elapsedTime;
	
	public static void main(String[] args) throws FileNotFoundException 
	{
		Scanner in;
		if (args.length ==  0) 
		{
		    in = new Scanner(System.in);
		}
		else
		{
			File fInputFile = new File(args[0]); 
			in = new Scanner(fInputFile);
		}
		//Read the directed input graph
		g = Graph.readGraph(in, true);
		iOriginalVerticesCount = g.verts.size();
		
		//Get MST using Edmond's Branching Algorithm
		timer();
		edmondBranchingMST();
		timer();
	}
	
	//Starts execution of Edmond's Branching Algorithm
	public static void edmondBranchingMST()
	{
		//Continue till all the active vertices in the graph are not reachable from root node by 0 weight edges
		boolean bAllVerticesReachable = areAllVerticesReachableFromRoot();
		while(bAllVerticesReachable == false)
		{
			//For each vertex find the incoming edge with the smallest weight
			for(int i=2; i< g.verts.size(); i++)
			{
				if(g.verts.get(i).activeVertex == true) //We don't want to process inactive vertices
				{
					reduceWeight(g.verts.get(i));				
				}							
			}
			bAllVerticesReachable = areAllVerticesReachableFromRoot();
		}
		
		if(iOriginalVerticesCount > 50)
		{
			System.out.println("Weight of MST:  " + MSTWeight);
		}
		else
		{
			//If any cycles were contracted, expand them now
			while(iContractedCycles > 0)
			{
				unContractCycles();
				iContractedCycles--;
			}
			
			System.out.println("Edges in MST:  ");
			for(int i=0; i< MST.size(); i++)
			{
				System.out.println(MST.get(i).toString());
			}
			System.out.println("Weight of MST:  " + MSTWeight);
		}		
	}	
	
	//Expands a cycle that was contracted before
	public static void unContractCycles()
	{
		ArrayList<Edge> MSTNew = new ArrayList<Edge>();
		for(int i=0; i< MST.size(); i++)
		{		
			if(MST.get(i).From.bContractedVertex == false && MST.get(i).To.bContractedVertex == false)
			{
				//Not a super vertex so need not be expanded
				MSTNew.add(MST.get(i));
			}
			else if(MST.get(i).From.bContractedVertex == false && MST.get(i).To.bContractedVertex == true)
			{
				//Edge from a vertex to a super vertex which needs to be expanded
				Vertex v = MST.get(i).From;
				Vertex v2= null;
				for(Edge e: v.Adj)
				{
					if(contractedVertices.get(iContractedCycles-1).contains(e.otherEnd(v)))
					{
						MSTNew.add(e);
						v2 = e.otherEnd(v);
						break;
					}
				}				
				//Take path from contracted cycle that starts from this vertex till we don't come back to it leave the last edge
				ArrayList<Edge> edgesWithoutCycle = rearrangeCycle(v2, contractedCycles.get(iContractedCycles-1));
				MSTNew.addAll(edgesWithoutCycle);
			}
			else if(MST.get(i).From.bContractedVertex == true && MST.get(i).To.bContractedVertex == false)
			{
				int iLength = contractedVertices.get(iContractedCycles-1).size();
				Vertex v = contractedVertices.get(iContractedCycles-1).get(iLength-1);
				for(Edge e3: v.Adj)
				{
					if(e3.otherEnd(v) == MST.get(i).To)
					{
						MSTNew.add(e3);
					}
				}
			}
		}
		
		MST = MSTNew;
	}
	
	//Rearranges the cycle edges so the first edge starts from input vertex
	//Also removes the last edge in the cycle so that the cycle is not formed again
	public static ArrayList<Edge> rearrangeCycle(Vertex v, ArrayList<Edge> cycle)
	{
		ArrayList<Edge> rearrangedCycle = new ArrayList<Edge>();
		rearrangedCycle = cycle;
		int iIndex = -1;
		for(int i=0; i< rearrangedCycle.size(); i++)
		{
			if(rearrangedCycle.get(i).From == v)
			{
				iIndex = i;
				break;
			}
		}
		ArrayList<Edge> rearrangedCycleTemp = new ArrayList<Edge>();
		for(int i=0; i < iIndex; i++)
		{
			rearrangedCycleTemp.add(rearrangedCycle.get(i));	
			rearrangedCycle.remove(i);
		}
		rearrangedCycle.addAll(rearrangedCycleTemp);
		rearrangedCycle.remove(rearrangedCycle.size()-1); //Remove last edge that forms the cycle
		return rearrangedCycle;
	}
	
	//Calculates the minimum  weight of all active incoming edges to the input vertex and then reduce the weight of all active
	//edges to the input vertex by that value
	public static void reduceWeight(Vertex v)
	{
		int iMin = 0;
		if(v.revAdj.isEmpty() && v.activeVertex == true)
		{
			iMin = 0; //If no incoming edges - set min weight to 0
			v.minWeight = iMin;
		}
		else
		{
			iMin = Integer.MAX_VALUE;
			if(v.activeVertex == true)
			{
				for(Edge e: v.revAdj)
				{
					if(e.activeEdge == true)
					{
						if(e.Weight < iMin)
						{
							iMin = e.Weight;
						}
					}				
				}				
				//Set the min value of the vertex as the minimum weight of all its incoming edges
				v.minWeight = iMin;
				MSTWeight = MSTWeight + iMin;
				//For each incoming edge to each vertex reduce the value of the weight by the smallest weight we just found
				for(Edge e: v.revAdj)
				{
					if(e.activeEdge == true)
					{
						e.Weight = e.Weight - v.minWeight;
					}
				}				
			}
		}
	}
	
	//Checks if all the active vertices are reachable from root vertex
	public static boolean areAllVerticesReachableFromRoot()
	{
		//Start the path from root vertex
		//Iterate over all vertices and find a zero weight edge connecting one vertex to other
		//If all vertices are part of the list then we can start expanding the vertices if they were contracted first
		//If the new vertex we found is already in the list - then we found a cycle in graph that needs to be contracted
		//If there is a vertex that was not reachable from root using 0 weight edge - we need to find a cycle including that vertex and contract it 
		Vertex root = g.verts.get(1);
		ArrayList<Vertex> vertices = new ArrayList<Vertex>(); //Contains vertices reachable from 1 by 0 edges
		ArrayList<Vertex> verticesToProcess = new ArrayList<Vertex>(); //Contains list of vertices yet to be processed for 0 weight edges possibilities
		Vertex vNext = root;
		verticesToProcess.add(root);
		while(!verticesToProcess.isEmpty())
		{
			vNext = verticesToProcess.remove(0);
			if(vNext.activeVertex == true)
			{
				for(Edge e: vNext.Adj)
				{
					if(e.Weight == 0 && e.activeEdge == true)
					{
						if(vertices.contains(e.otherEnd(vNext)))
						{
							//cycle found
							break;
						}
						else
						{
							//Not a cycle yet so continue to find next zero weight edge from the new vertex we added to the list
							vertices.add(e.otherEnd(vNext)); //Add the vertex to the list of vertices reachable from root using 0 weight edges
							verticesToProcess.add(e.otherEnd(vNext)); //Add vertices to list of nodes that still need to be processed 
																	  //for 0 weight edges to other vertices in graph
						}						
					}
				}
			}
		}	
		int iActiveVertices = countActiveVertices();
		if(iActiveVertices-1 == vertices.size()) // subtract 1 since root is not included in the vertices list
		{
			//All active vertices are reachable from vertex 1 (root)
			findMSTpath(vertices);
			return true; //Start expanding vertices if they were contracted before
		}
		else
		{
			ArrayList<Vertex> verticesToEvaluateForCycles = new ArrayList<Vertex>();
			//Try to contract a zero weight cycle from a randomly selected vertex
			//Continue till no zero weight cycle is found and then contract it to a super vertex
			verticesToEvaluateForCycles = findRemainingVertices(vertices);
			findZeroWeightCycle(verticesToEvaluateForCycles);
			return false;//Since all vertices are not reachable from root by 0 weight edges return false
		}		
	}
	
	//Finds the MST path of zero weight edges from root to all the vertices in the input list 
	public static void findMSTpath(ArrayList<Vertex> vertices)
	{		
		Vertex root = g.verts.get(1);
		Vertex vNext = root;
		Edge e = null;
		for(int i=0; i< vertices.size(); i++)
		{
			e = returnZeroWeightEdgeOutgoing(vNext);
			if(e.otherEnd(vNext) == vertices.get(i))
			{
				MST.add(e);
				vNext = e.otherEnd(vNext);
			}
		}
	}
	
	//Finds a zero weight cycle between the vertices passed as argument
	//If a cycle is found then the cycle is contracted to a super vertex
	public static void findZeroWeightCycle(ArrayList<Vertex> verticesList)
	{
		ArrayList<Edge> zeroWtCycle = new ArrayList<Edge>();
		ArrayList<Vertex> zeroVertices = new ArrayList<Vertex>();		
		Vertex startingVertex = verticesList.get(0);
		zeroVertices.add(startingVertex);
		Edge nextEdge = returnZeroWeightEdge(startingVertex);
		Vertex nextVertex = startingVertex;
		while(nextEdge != null)
		{
			zeroWtCycle.add(nextEdge);
			nextVertex = nextEdge.otherEnd(nextVertex);			
			if(zeroVertices.contains(nextVertex))
			{
				//We found a zero weight cycle
				contractedCycles.add(zeroWtCycle); //Saving this cycle as we need to expand it later
				ArrayList<Edge> zeroWtCycleNew = removeEdgesNotInCycle(zeroWtCycle);
				iContractedCycles++;
				contractVertex(nextVertex,zeroWtCycleNew);
				break;
			}
			else
			{
				zeroVertices.add(nextVertex);
				nextEdge = returnZeroWeightEdge(nextVertex);
			}
		}//End of while(!nextEdge.equals(null))
	}
	
	//Removes the edges that are start from one of the vertices in cycle to a vertex not part of cycle
	//Also removes the edges that lead to one of the vertices in cycle from a vertex not part of cycle
	public static ArrayList<Edge> removeEdgesNotInCycle(ArrayList<Edge> cycle)
	{
		ArrayList<Edge> edges = new ArrayList<Edge>();
		edges = cycle;
		int[][] arr = new int[cycle.size()][2];
		int iVertex1 = 0;
		int iVertex2 = 0;
		boolean bNew1 = true;
		boolean bNew2 = true;
		int k=-1;
		for(int i=0; i< cycle.size(); i++)
		{
			iVertex1 = cycle.get(i).From.name;
			iVertex2 = cycle.get(i).To.name;
			bNew1 = true;
			bNew2 = true;
			for(int j=0; j< arr.length; j++)
			{
				if(arr[j][0] == iVertex1)
				{
					bNew1 = false;
					arr[j][1]++;
				}
				
				if(arr[j][0] == iVertex2)
				{
					bNew2 = false;
					arr[j][1]++;
				}
			}
			if(bNew1 == true)
			{
				k++;
				arr[k][0] = iVertex1;
				arr[k][1] = 1;
			}
			
			if(bNew2 == true)
			{
				k++;
				arr[k][0] = iVertex2;
				arr[k][1] = 1;
			}
		}
		
		Edge e;
		for(int i=0; i< arr.length; i++)
		{
			if(arr[i][1] == 1) //this vertex is not part of cycle - remove that edge
			{
				for(int j=0; j< edges.size(); j++)
				{
					e = edges.get(j);
					if(e.From.name == arr[i][0] || e.To.name == arr[i][0])
					{
						edges.remove(j);
					}
				}
			}
		}
		return edges;
	}

	//Contract the input cycle given
	public static void contractVertex(Vertex vStartVertex, ArrayList<Edge> cycle)
	{
		//To contract the cycle - create a new vertex in graph g
		Vertex v = new Vertex(g.verts.size());
		g.verts.add(v);
		v.activeVertex = false;
		v.bContractedVertex = true;
		
		//Mark the edges in cycle as inactive
		for(Edge e: cycle)
		{
			e.activeEdge = false;
		}
		
		//Get the list of all vertices part of the zero weight cycle
		ArrayList<Vertex> zeroVertices = new ArrayList<Vertex>();
		for(Edge e: cycle)
		{
			zeroVertices.add(e.From);
		}
		//For each zero vertex if the incoming edge comes from another zero vertex
		//Mark that edge inactive and add it to the list of contracted edges,  edge can be a non zero weight edge or otherwise
		//This edge is not be part of cycle
		//Now go through all the edges for each active vertex in  the graph
		//If an edge is from one non 0 vertex to another non 0 vertex - keep it active - nothing to do here
		//If an edge is from one non 0 vertex to a 0 vertex - add that edge to new vertex and make the original edge inactive 
			//use updated weights not the original weight of the edge
		//If the edge is from 0 vertex to non 0 vertex the add that edge to non 0 vertex coming from new vertex 
			//make the original edge inactive and the new edge active - use updated weight not the original weight of the edge
		Vertex u;
		ArrayList<Edge> updatedEdges = new ArrayList<Edge>();
		for(int i=1; i< g.verts.size(); i++)
		{
			u = g.verts.get(i);
			if(u.activeVertex == true)
			{
				if(zeroVertices.contains(u))
				{
					for(Edge e: u.Adj)
					{
						//This is an edge from one zero vertex to another - mark this as inactive edge
						if(zeroVertices.contains(e.otherEnd(u)))
						{
							e.activeEdge = false;
						}
						else
						{
							//The edge is from one zero vertex to a non zero vertex - we need to add this edge to the new vertex
							Edge newEdge = new Edge(v, e.To, e.Weight);
							v.Adj.add(newEdge);
							e.To.revAdj.add(newEdge);
							e.activeEdge = false;
						}
					}
					for(Edge e: u.revAdj)
					{
						//This is an edge from one zero vertex to another - mark this as inactive edge
						if(zeroVertices.contains(e.otherEnd(u)))
						{
							e.activeEdge = false;
						}
						else
						{
							//The edge is from one zero vertex to a non zero vertex - we need to add this edge to the new vertex
							Edge newEdge = new Edge(e.From, v, e.Weight);
							e.From.Adj.add(newEdge);
							v.revAdj.add(newEdge);
							e.activeEdge = false;
						}
					}
				}
				else //If not a zero vertex then check if any edge leads to a zero vertex
				{
					for(Edge e1: u.Adj)
					{
						if(zeroVertices.contains(e1.otherEnd(u)))
						{
							Edge newEdge = new Edge(u, v, e1.Weight);
							newEdge.activeEdge = true;
							updatedEdges.add(newEdge);
							v.revAdj.add(newEdge);
							e1.activeEdge = false;
						}
						else //If edge is from non zero vertex to non zero vertex - do nothing
						{
							continue;
						}
					}
					u.Adj.addAll(updatedEdges);
					updatedEdges.clear();
					for(Edge e1: u.revAdj)
					{
						if(zeroVertices.contains(e1.otherEnd(u)))
						{
							Edge newEdge = new Edge(v, u, e1.Weight);
							newEdge.activeEdge = true;
							v.Adj.add(newEdge);
							updatedEdges.add(newEdge);
							e1.activeEdge = false;
						}
						else //If edge is from non zero vertex to non zero vertex - do nothing
						{
							continue;
						}
					}
					u.revAdj.addAll(updatedEdges);
					updatedEdges.clear();
				}
			}
		}
		//Make the new vertex active
		v.activeVertex = true;
		//Make all zero vertices inactive
		contractedVertices.add(zeroVertices);
		for(int i=0; i< zeroVertices.size(); i++)
		{
			zeroVertices.get(i).activeVertex = false;
		}
	}
	
	//Returns the list of vertices not reachable from root vertex by 0 weight edges
	public static ArrayList<Vertex> findRemainingVertices(ArrayList<Vertex> verticesList)
	{
		ArrayList<Vertex> verticesRemaining = new ArrayList<Vertex>();
		for(int i=2; i< g.verts.size(); i++)
		{
			if(!verticesList.contains(g.verts.get(i)) && g.verts.get(i).activeVertex == true)
			{
				verticesRemaining.add(g.verts.get(i));
			}
		}		
		return verticesRemaining;
	}
	
	//Returns a zero weight outgoing edge for the input vertex
	//If no such edge exists then it returns null
	public static Edge returnZeroWeightEdgeOutgoing(Vertex v)
	{
		if(v.Adj.size() > 0)
		{
			for(Edge e: v.Adj)
			{
				if(e.Weight == 0 && e.activeEdge == true)
				{
					return e;
				}				
			}
		}	
		return null;
	}
	
	//Finds a zero weight incoming edge to the vertex v and returns that edge
	//If no zero weight incoming edge or no incoming edges at all to the vertex - returns null
	public static Edge returnZeroWeightEdge(Vertex v)
	{
		if(v.revAdj.size() > 0)
		{
			for(Edge e: v.revAdj)
			{
				if(e.Weight == 0 && e.activeEdge == true)
				{
					return e;
				}				
			}
		}	
		return null;
	}
	
	//Counts and returns the number of active vertices in the graph
	public static int countActiveVertices()
	{
		int iCountActiveVertices = 0;
		for(int i=1;i< g.verts.size(); i++)
		{
			if(g.verts.get(i).activeVertex == true)
			{
				iCountActiveVertices++;
			}
		}
		return iCountActiveVertices;
	}

	public static void timer() 
	{
		if (phase == 0) 
		{
			startTime = System.currentTimeMillis();
			phase = 1;
		} 
		else 
		{
			endTime = System.currentTimeMillis();
			elapsedTime = endTime - startTime;
			System.out.println("Time: " + elapsedTime + " msec.");
			memory();
			phase = 0;
		}
	}

	public static void memory() 
	{
		long memAvailable = Runtime.getRuntime().totalMemory();
		long memUsed = memAvailable - Runtime.getRuntime().freeMemory();
		System.out.println("Memory: " + memUsed / 1000000 + " MB / " + memAvailable / 1000000 + " MB.");
	}

}


