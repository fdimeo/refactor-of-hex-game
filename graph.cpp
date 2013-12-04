//====================================================================
//
//  THIS PROGRAM MUST BE COMPILED WIH THE C++11 extensions
//
//  Compiled using Gnu compiler version 4.8
//  
//  g++ graph.cpp -std=C++11 -o graph
//
//====================================================================
#include <iostream>
#include <map>
#include <vector>
#include <list>
#include <array>
#include <ctime>    // standard C library
#include <cstdlib>  // standard C library
#include <string>
#include <chrono>
#include <random>
#include <algorithm>

#include <fstream>

// forward class declarations
class graphPoint;
class Graph;
class HexPlayer;
class monteCarlo;

// node color class
enum class nodecolor{ 
   NONE, RED, BLUE };

// the playertype definition
enum class playertype { 
   HUMAN, COMPUTER };

typedef nodecolor playercolor;    // the node colors are also the player colors
typedef unsigned int movevalue;   // the resultant value of a generic move (range is 0-worst to 100-best)

typedef std::vector< unsigned int > uintvec;

// output operator for the playercolor
std::ostream &operator<<(std::ostream &out, playercolor color)
{
   out << ( color == playercolor::RED ? "\e[1;31mRED\e[0m" : "\e[1;36mBLUE\e[0m");
   return out;
}



//
// the hexGame class definition
//
class hexGame 
{

 private:
   unsigned int m_boardsize;        // the size of the hex playing board
   Graph *m_pBoard;                 // a pointer to an instance of the board
   bool mstIncludesColumn(unsigned int column, std::vector< int > *pMstVector);
   bool mstIncludesRow(unsigned int row, std::vector< int > *pMstVector);
   void construct_board();

 public:
   hexGame();
   ~hexGame();
   hexGame(unsigned int gamesize);
   unsigned int getGameSize();
   bool isMoveLegal(unsigned int row, unsigned int col);               // returns true if move is legal
   bool makeGameMove(HexPlayer &p, unsigned int row, unsigned int col);    // returns true if move was legal
   bool makeGameMove(HexPlayer &p, unsigned int node);
   int getTileNumber(unsigned int row, unsigned int col); // compute the absolute position number given the row and node number in that row
   int getTileRow(unsigned int nodeNumber);
   int getTileColumn(unsigned int nodeNumber);
   nodecolor getTileColor(unsigned int nodeNumber);
   bool getAllFreeTiles(std::vector< unsigned int > &pFreeTiles);
   void forceTileColor(unsigned int node, playercolor color);

   // checks to see if the player specified has won
   bool checkForWinner(HexPlayer *player);

   friend std::ostream& operator<<(std::ostream &out, hexGame &game);

};



// A generic player.  Type of game unknown.
class Player 
{
 protected:
   unsigned int index;                              // the assigned index for the player (0-none)
   std::string m_name;                              // the name of the player

 public:
   Player();
   Player(std::string name);
   virtual void makeMove() = 0;                     // pure virtual function for making a move
   unsigned int getPlayerIndex();
   void setPlayerIndex(unsigned int player_index);
   
};

// A "game of hex" player
//
class HexPlayer:public Player 
{
protected:
   playercolor m_color;                            // playercolor::RED or playercolor::BLUE
   hexGame &m_game;                                // a reference to the gameboard being used in this game
   playertype m_type;                              // the type of player (human, computer)
   unsigned int m_lastMove;                        // the node number (successfully) chosen last

public:
   HexPlayer(hexGame &game, std::string name);
   playercolor getPlayerColor();
   void setPlayerColor(playercolor color);
   virtual void makeMove() = 0;                    // pure virtual move function
   void setPlayerType(playertype type);
   playertype getPlayerType();
   void setLastMove(unsigned int nodenumber);
   unsigned int getLastMove();
};

// A computer hex game player
class ComputerHexPlayer:public HexPlayer
{
private:
   monteCarlo *mc_sim;                         // and instance of the monte carlo simulator

public:
   ComputerHexPlayer(hexGame &game );
   void makeMove();                             // makes a move appropriate for a computer player
   
};

// a human hex game player
class HumanHexPlayer:public HexPlayer
{
public:
   HumanHexPlayer(hexGame &game);
   void makeMove();                            // makes a move appropriate for a human player
   
};

std::ostream &operator<<(std::ostream &out, std::vector< unsigned int > &vec_of_ints)
{
   out << "ints in vector are:" << std::endl;

   for(std::vector< unsigned int >::iterator it=vec_of_ints.begin(); it<vec_of_ints.end(); it++)
   {
      out << *it << std::endl;
   }

   return out;

}
//
// monte carlo class definition
//
class monteCarlo
{
private:
   static const int NUM_OF_TRIALS=1000;
   hexGame *m_p_mcsim_game;                            // our private (hex board) to run the simulation on
   hexGame &m_hexgame;                                 // a reference to the hex game that we're running the simulation against
   std::vector<double> trials_results;                 // the percentage win, one for each possible move on the board
   std::vector< unsigned int > next_moves;             // a vector of next moves for evaluation
   std::vector< unsigned int > all_possible_moves;     // one move for each possible move on the board
   
public:
   monteCarlo(hexGame &m_hexgame);                    // a reference to the hex board that we're running the mc simulation on

   // A functor because...well...why not? :)  
   //  Runs the simulation logic and returns the best node (zero based) for the next move for the game specified
   unsigned int operator()(HexPlayer &p);        

};


monteCarlo::monteCarlo(hexGame &hexgame):m_hexgame(hexgame)
{
   int gamesize = m_hexgame.getGameSize();
   std::cout << "Creating a simulation hex board of dimension " << gamesize << std::endl;
   m_p_mcsim_game = new hexGame(gamesize);   // a new simulation game
   std::cout << "simulation hex board of dimension " << m_p_mcsim_game->getGameSize() << " created " << std::endl;
   
}

//===========================================================================================
//  Monte carlo simulaton logic
//-------------------------------------------------------------------------------------------
//
// copy the board into our class
//
// find all the empty nodes in the player game, list them in all_possible_moves
// for each node in all_possible_moves
//    place that node on the board AND put it in next_moves location
//    for each iteration in NUM_OF_TRIALS
//       shuffle the vector copy of next possible moves
//       place each tile in the vector copy of next moves on the board, alternate red and blue
//       see who won
//       compute the trial results and keep the average in trails_results
// 
//  look through each of the trials_results for the largest number
//  get the associated random_next_move and return it as the answer
//
#define MONTE_CARLO_DEBUGGED

unsigned int monteCarlo::operator()(HexPlayer &p)        
{
   unsigned int computer_move;

#ifdef MONTE_CARLO_DEBUGGED

   std::cout << "Running MC simulation using simulated game of size " << m_p_mcsim_game->getGameSize() << std::endl;

   std::cout <<  m_hexgame << std::endl;

   unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

   // copy the graph to our simluation graph
   for(unsigned int node=0; node < (m_p_mcsim_game->getGameSize() * m_p_mcsim_game->getGameSize()); node++)
   {
      std::cout << "setting node " << node << " to " << m_hexgame.getTileColor(node) << std::endl;
      m_p_mcsim_game->forceTileColor(node, m_hexgame.getTileColor(node));
   }

   // find all the moves we could make
   m_p_mcsim_game->getAllFreeTiles(all_possible_moves);

   // print them out
   //   std::cout << all_possible_moves << std::endl;

   // create a vector of the remaining moves 
   uintvec remaining_moves = all_possible_moves;

   // try each move as the next move 
   for(uintvec::iterator it=all_possible_moves.begin(); it!=all_possible_moves.end(); it++)
   {
      unsigned int next_move = *it;

      // remove the move from the possible moves
      remaining_moves.erase(remaining_moves.begin());

      // print them out
      //      std::cout << remaining_moves << std::endl;

      m_p_mcsim_game->forceTileColor(next_move, p.getPlayerColor());
      
      next_moves.push_back(next_move);  // this is the next move we are trialing

      nodecolor tile_color = ((p.getPlayerColor() == playercolor::RED) ? playercolor::BLUE : playercolor::RED);

      // do this simulation NUM_OF_TRIALS times, compute the win probability as we go
      for(int num_trials=0; num_trials<1; num_trials++)
      {
         // shuffle the remaining moves
         std::shuffle (remaining_moves.begin(), remaining_moves.end(), std::default_random_engine(seed));
         
         // now fill the board with alternating colors
         for (uintvec::iterator rit=remaining_moves.begin(); rit!=remaining_moves.end(); rit++)
         {
            m_p_mcsim_game->forceTileColor(*it, tile_color);
            tile_color = ((tile_color == playercolor::RED) ? playercolor::BLUE : playercolor::RED);
         }

         std::cout << "==============================" << std::endl;
         std::cout << *m_p_mcsim_game << std::endl;
         std::cout << "==============================" << std::endl;
         
         char dbg;
         std::cout << "hit c to continue";
         std::cin >> dbg;

         // is this a winner for the player?
         

      }
   }

   return 0;  // for now

#else

   while(true)
   {
      // a number between 0 and the board size
      computer_move = rand() % (m_hexgame.getGameSize() * m_hexgame.getGameSize());

      std::cout << "computer_move is node number " << computer_move << std::endl;

      if (m_hexgame.isMoveLegal(m_hexgame.getTileRow(computer_move), m_hexgame.getTileColumn(computer_move)))
      {
         // we've picked our move
         break;
      }
   }

#endif

   std::cout << p.getPlayerColor() << " chooses tile " << m_hexgame.getTileColumn(computer_move)+1 << ", " << m_hexgame.getTileRow(computer_move)+1 << std::endl;

   // return the move that we've chosen
   return computer_move;


}



class ShortestPathAlgo
{
private:
   std::list<unsigned int> *pathList;
   int pathCost;  // the path cost, or (-1) if no path exists

public:

   ShortestPathAlgo();
   ~ShortestPathAlgo();

   // returns a count of the nodes
   unsigned int verticies(Graph &G);

   // returns the cost of the path (or -1 if no path exists)
   int path_size( Graph &G, unsigned int originNode, unsigned int destNode );

   // returns a std::list pointer with the path
   std::list<unsigned int> *path( Graph &G, unsigned int originNode, unsigned int destNode);

   // this helps print the path list
   friend std::ostream &operator<< (std::ostream &cout, std::list<unsigned int> *path);


};
//-------------------------------------------------------------------------------------------------------
//  A class defining an entire graph, which is comprised of graphPoints with edges to other graphPoints
//-------------------------------------------------------------------------------------------------------
//
class Graph
{

private:
   std::map< int, graphPoint* > graphNodes;// a map of all graphPoints (i.e. nodes, vertices) in the graph
   unsigned int m_totalNumVerticies;       // the total number of vertices (nodes) in this graph
   unsigned int m_totalNumEdges;           // the total number of edges in this graph
   int m_originNode;                       // the node number of the origin (-1 if not assigned)

public:
   Graph();                                // constructor for normal graph
   Graph(unsigned int hexGraphDimension);  // constructor for a hex graph of the given dimension
   Graph( char* fileName);                 // constructor for a graph using edges defined in a local file
   void addNode(unsigned int nodeNumber);
   void removeNode(unsigned int nodeNumber);
   void setNodeValue(unsigned int nodeNumber, unsigned int cost);
   int getNodeValue(unsigned int nodeNumber);
   void setNodeColor(unsigned int nodeNumber, nodecolor color);
   nodecolor getNodeColor(unsigned int nodeNumber);
   bool isNodeInGraph(unsigned int nodeNumber);
   void addEdge(unsigned int sourceNodeNumber, unsigned int destNodeNumber, unsigned int edgeWeight);
   bool hasEdge(unsigned int sourceNodeNumber, unsigned int destNodeNumber);

   void deleteEdge(unsigned int sourceNodeNumber, unsigned int destNodeNumber);
   void makeOriginNode(unsigned int nodeNumber);  // origin nodes have a totalcost of 0 and have been visited
   int modifyEdge(unsigned int destNodeNumber, unsigned int weight);
   int getEdgeValue(unsigned int sourceNodeNumber, unsigned int destNodeNumber);
   int setEdgeValue(unsigned int sourceNodeNumber,unsigned int destNodeNumber, unsigned int weight );
   unsigned int getNodeCount(void);
   unsigned int getEdgeCount(void);
   bool isNodeVisited(unsigned int nodeNumber);
   void setNodeVisited(unsigned int nodeNumber);
   void resetNodeVisited(unsigned int nodeNumber);
   void resetNodeVisitedAll();


   void doDijkstra( unsigned int originNode, unsigned int destNode, std::list<unsigned int> *pathResult, int &pathCost);
   bool doPrim( unsigned int originNode, nodecolor color,  std::vector< int > *pVectorMSTNodeResult, bool printSolution);

   void printGraph();
};


//-------------------------------------------------------------------------------------------------------
//  A class defining a graphPoint, which may become included in a graph
//-------------------------------------------------------------------------------------------------------
//
class graphPoint 
{

   friend class Graph;

private:
   unsigned int m_nodeNumber;                    // a unique identifier for this node
   std::map<unsigned int, unsigned int> m_edges; // a vector of all edges from the node (other node number, cost)
   unsigned int m_viaNode;                       // the "from node" to this node for the m_totalCost recorded
   int    m_totalCost;                           // total path cost for this instance
   bool   m_visited;                             // has this instance been visited in the algorthim?
   unsigned int m_numEdges;                      // this number of edges for this instance
   unsigned int m_hexRowNumber;                  // the hex row number ( if this is a hex graph, or (-1) if not )
   nodecolor m_nodeColor;                        // NONE, RED, or BLUE

public:
   graphPoint( unsigned int nodeNumber, int cost, bool visited);
   void printGraphPoint();
   void createEdge(unsigned int dest_node, unsigned int weight);
   int deleteEdge(unsigned int dest_node);
   int setEdgeValue(unsigned int sourceNodeNumber,unsigned int NodeNumber, unsigned int weight );
   int getEdgeValue(unsigned int sourceNodeNumber );
   int modifyEdge(unsigned int dest_node, unsigned int weight);
   int getPointCost ();
   void setPointCost (int cost);
   bool modifyPoint (int cost, bool visited);
   void setVisited ();
   void resetVisited ();
   bool getVisited ();
   void cleanNode();


};

 
//*****************************************************************
//**
//** graphPoint methods
//**
//*****************************************************************

graphPoint::graphPoint( unsigned int nodeNumber, int cost = (-1), bool visited = false )
{
   m_totalCost = cost;          // the accumulated cost of this node
   m_visited = visited;         // not visited when created (but can be overriden for the source node)
   
   m_nodeNumber = nodeNumber;   // this node's number
   m_viaNode = nodeNumber;      // the node number that this node was reached from (start with self)
   m_numEdges = 0;              // there are no edges to start
   m_hexRowNumber = -1;         // this is not a hex node element
   m_nodeColor = nodecolor::NONE;

}


bool graphPoint::modifyPoint (int cost, bool visited = false)
{
   m_totalCost = cost;
   m_visited = visited;
}

void graphPoint::setVisited ()
{
   m_visited = true;
}

void graphPoint::resetVisited ()
{
   m_visited = false;
}


bool graphPoint::getVisited ()
{
   return m_visited;
}


int graphPoint::getPointCost ()
{
   return m_totalCost;
}

void graphPoint::setPointCost (int cost)
{
   m_totalCost = cost;
}

void graphPoint::cleanNode()
{
   m_viaNode = 0;
   m_totalCost = -1;
   m_visited = false;
}


void graphPoint::printGraphPoint()
{
   std::cout << "Graph point #" << m_nodeNumber
             << ((m_totalCost == 0) ? " (ORIGIN)" : "")
             << " has a cost of " 
             << m_totalCost 
             << " and has" << (m_visited ? "" : " not") << " been visited" << std::endl;
   
   if(m_numEdges) std::cout << "Edge at:" << std::endl;
   
   for(std::map<unsigned int, unsigned int>::iterator it=m_edges.begin(); it != m_edges.end(); ++it)
   {
      std::cout << "-- to node:" << it->first << " (" << it->second << ")" << std::endl;
   } 
}

// create a new edge to "dest_node" with a cost of "weight"
void graphPoint::createEdge(unsigned int dest_node, unsigned int weight)
{
      std::map<unsigned int, unsigned int>::iterator it;

      // replace any duplicate entry with a new weight
      if((it = m_edges.find(dest_node)) != m_edges.end())
      {
         // std::cout << "duplicate edge \"dest node\"found: "
         //           << it->first << ", replacing " << it->second 
         //           << " with " << weight << std::endl;

         it->second = weight;
         return;
      }

      // otherwise, insert a new edge in the set of all edges from this node
      else
      {
         m_edges.insert(std::pair<unsigned int,unsigned int>(dest_node, weight));
         m_numEdges++;
     }
}

// remove the edge to "dest_node"
int graphPoint::deleteEdge(unsigned int dest_node)
{
   int retval = -1;
   std::map<unsigned int, unsigned int>::iterator it;
   
   // find and delete the edge.  If not found, do nothing
   if((it = m_edges.find(dest_node)) != m_edges.end())
   {
      std::cout << "\"dest node\"found: "
                << it->first << ", erasing edge to it " << it->second 
                << std::endl;
      
      m_edges.erase(it);
      m_numEdges--;
      retval = 0;
   }

   return retval;
   
}

// modify the edge to "dest_node"
//
// return -1 if error, 0 if ok.
//
int graphPoint::modifyEdge(unsigned int dest_node, unsigned int weight)
{
   int retval = -1;
   std::map<unsigned int, unsigned int>::iterator it;
   
   // find and delete the edge.  If not found, do nothing
   if((it = m_edges.find(dest_node)) != m_edges.end())
   {
      // std::cout << "\"dest node\"found: "
      //           << it->first << ", modifying the edge to " << it->second 
      //           << std::endl;
      
      it->second = weight;
      return 0;
   }

   return -1;

}

// get the edge to "dest_node"
//
// return -1 if error, edge weight if ok
//
int graphPoint::getEdgeValue(unsigned int node)
{
      int retval = -1;
      std::map<unsigned int, unsigned int>::iterator it;

      // std::cout << "Finding edge cost from " << m_nodeNumber << " to " << node << std::endl;

      // find the edge.  If not found, return a cost of -1 (this mean infinity)
      if((it = m_edges.find(node)) != m_edges.end())
      {
         // std::cout << "Found a cost of : " << it->second << std::endl;

         retval = (it->second);
      }

      // std::cout << "done" << std::endl;

      return retval;

}



//*****************************************************************
//**
//** Graph methods
//**
//*****************************************************************
//
// the constructor for "normal", free form graphs
Graph::Graph()
{
   m_totalNumVerticies = 0;
   m_totalNumEdges = 0;
}

Graph::Graph(char *fileName)
{

   m_totalNumVerticies = 0;
   m_totalNumEdges = 0;


   int graphSize;
   
   // read the graph from a file.
   // the first int is the size, and followed by tuples of int node1, int node2, int cost
   // open the graph file for read only
   std::fstream filestream(fileName, std::fstream::in );

   filestream >> graphSize;

   std::cout << " graph size is " << graphSize << " nodes" << std::endl;

   while(true)
   {    
      unsigned int node1, node2, cost;
      
      // collect a edge data
      filestream >> node1 >> node2 >> cost;

      // only read with the input is valid
      if(filestream.eof()) break;

      std::cout << "node1: " << node1 << " node2: " << node2 << " cost: " << cost << std::endl;
      addNode(node1);       // these addNode calls are idempotent 
      addNode(node2);
      addEdge(node1, node2 , cost);
   }
}


// add a node to the graph
void Graph::addNode(unsigned int nodeNumber)
{
   std::map<int, graphPoint* >::iterator it = graphNodes.find(nodeNumber);

   if(it == graphNodes.end())
   {
      graphNodes[nodeNumber] = new graphPoint(nodeNumber);
      m_totalNumVerticies++;
      graphNodes[nodeNumber]->resetVisited();
      m_originNode = -1;
   }
}

void Graph::addEdge(unsigned int sourceNodeNumber, unsigned int destNodeNumber, unsigned int edgeWeight)
{
   std::map<int, graphPoint* >::iterator it_s = graphNodes.find(sourceNodeNumber);
   std::map<int, graphPoint* >::iterator it_d = graphNodes.find(destNodeNumber);

   if( it_s != graphNodes.end())
   {
      it_s->second->createEdge(destNodeNumber, edgeWeight);
      m_totalNumEdges++;
   }
}

bool Graph::hasEdge(unsigned int sourceNodeNumber, unsigned int destNodeNumber)
{
    return (graphNodes.find(sourceNodeNumber)->second->getEdgeValue(destNodeNumber) > 0);
}


//  make the specified node the origin
void Graph::makeOriginNode(unsigned int nodeNumber)
{
   std::map<int, graphPoint* >::iterator it = graphNodes.find(nodeNumber);

   // To be complete, we should check if m_originNode is -1 before we set the origin.
   // If m_originNode != -1, then the user had already set the origin, and is modifying it.
   // In that case, we should remove "orgin-ness" from the originally set origin point
   // For now, assume that the origin will not be changed once it's set.

   if( it != graphNodes.end())
   {
      it->second->modifyPoint (0, true);
      m_originNode = nodeNumber;
   }
}

void Graph::setNodeValue(unsigned int nodeNumber, unsigned int cost)
{
   std::map<int, graphPoint* >::iterator it = graphNodes.find(nodeNumber);

   if( it != graphNodes.end())
   {
      it->second->setPointCost (cost);
   }
}

int Graph::getNodeValue(unsigned int nodeNumber)
 {
    int retval = -1;
    std::map<int, graphPoint* >::iterator it = graphNodes.find(nodeNumber);

    if( it != graphNodes.end())   
    {
      retval = it->second->getPointCost ();
    }

    return retval;
 }

bool Graph::isNodeInGraph(unsigned int nodeNumber)
{
   
   int retval = false;
   std::map<int, graphPoint* >::iterator it = graphNodes.find(nodeNumber);
   
   if( it != graphNodes.end())   
   {
      retval = true;
   }
   
   return retval;
}


//returns -1 if not found
int Graph::getEdgeValue(unsigned int sourceNodeNumber,unsigned int destNodeNumber)
{
   unsigned int retval = -1;
   std::map<int, graphPoint* >::iterator it = graphNodes.find(sourceNodeNumber);

   if( it != graphNodes.end())
   {
      retval = it->second->getEdgeValue (destNodeNumber);
   }

   return retval;
}

bool Graph::isNodeVisited(unsigned int nodeNumber)
{

   bool retval = false;
   std::map<int, graphPoint* >::iterator it = graphNodes.find(nodeNumber);

   if( it != graphNodes.end())
   {
      retval = it->second->getVisited();
   }

   return retval;

}

void Graph::setNodeVisited(unsigned int nodeNumber)
{
   std::map<int, graphPoint* >::iterator it = graphNodes.find(nodeNumber);

   if( it != graphNodes.end())
   {
      it->second->setVisited();
   }

}

void Graph::resetNodeVisited(unsigned int nodeNumber)
{
   std::map<int, graphPoint* >::iterator it = graphNodes.find(nodeNumber);

   if( it != graphNodes.end())
   {
      it->second->resetVisited();
   }

}


void Graph::resetNodeVisitedAll()
{
   for(std::map<int, graphPoint* >::iterator it = graphNodes.begin(); it != graphNodes.end(); ++it)
   {
      it->second->resetVisited();
   }
}




enum nodecolor Graph::getNodeColor(unsigned int nodeNumber)
{
   std::map<int, graphPoint* >::iterator it = graphNodes.find(nodeNumber);

   if( it != graphNodes.end())
   {
      return it->second->m_nodeColor;
   }

   return nodecolor::NONE;

}

void Graph::setNodeColor(unsigned int nodeNumber, enum nodecolor color)
{
   std::map<int, graphPoint* >::iterator it = graphNodes.find(nodeNumber);

   if( it != graphNodes.end())
   {
      // only allow a color change if the node is already no color, or we are trying to
      //  reset the color of the tile back to NONE.  Once a tile is a play color, it
      //  can't be flipped back.  This allows us to call this method even when a tile
      //  has been taken by a player.
      if((it->second->m_nodeColor == nodecolor::NONE) || (color == nodecolor::NONE))
      {
         it->second->m_nodeColor = color;
      }
   }
}

unsigned int Graph::getNodeCount()
{
   return m_totalNumVerticies;
}

unsigned int Graph::getEdgeCount()
{
   return m_totalNumEdges;
}

void Graph::printGraph()
{
   unsigned int retval = -1;
   

   for(std::map<int, graphPoint* >::iterator it = graphNodes.begin(); it != graphNodes.end(); ++it)
   {
      it->second->printGraphPoint();
   }

   std::cout << "TOTAL NODES: " << getNodeCount() << "\tTOTAL EDGES: " << getEdgeCount() << "\n" << std::endl;
 

}


// a type which defines a conceptual two dimensional vector of ints (or at least can be addressed as such)
typedef typename std::vector< std::vector<int> > mst_result;

typedef typename std::map<unsigned int, unsigned int> edge_type;

// a type which describes a node in the graph
typedef typename std::map< int, graphPoint* > node_type;

const int SRC_NODENUM_IDX = 0;
const int DST_NODENUM_IDX = 1;
const int EDGEWEIGHT_IDX = 2;


//============================================================================================
//
// Graph member method: doPrim
// arguments: The origin node to start the MST computation from (unsigned int originNode)
// returns: Nothing
//
// Computed the Minimum Spanning tree solution for a Graph already instantiated in this class.
// The result of executing this method is to compute and print out the MST solution
//
//============================================================================================
bool Graph::doPrim( unsigned int originNode, nodecolor color = nodecolor::NONE, std::vector< int > *pSolutionVector = NULL, bool printSolution = true)
{

   int i;
   bool retval = false;

   if(getNodeColor(originNode) != color) return false;

   // initialize the 2dimenstional vector of [numNodes][2] ints. (nodenumber, weight) will
   //  be stored as we compute the prim solution
   mst_result mst_for_graph(m_totalNumVerticies, std::vector<int> (3));

   // initialize the solution 
   for(i=0; i<m_totalNumVerticies; i++)
   {
      mst_for_graph[i][SRC_NODENUM_IDX]=(-1);
      mst_for_graph[i][DST_NODENUM_IDX]=(-1);
   }
   
   // clear the visited indicators in all nodes of the graph
   resetNodeVisitedAll();

   // add the origin node to the solved set
   mst_for_graph[0][SRC_NODENUM_IDX] = originNode;
   mst_for_graph[0][DST_NODENUM_IDX] = originNode;
   mst_for_graph[0][EDGEWEIGHT_IDX] = 0;
   if(pSolutionVector) pSolutionVector->push_back(mst_for_graph[0][DST_NODENUM_IDX]);

   setNodeVisited(originNode);


   int solution_points_found = 1;  // the number of nodes in the MST solution thus far


    // now build the MST while iterating through the graph
   for(i=(m_totalNumVerticies-1); i>0; i--)
   {
      unsigned int lowest_cost_edge_this_iteration = 100; // the lowest cost and node for this iteration
      unsigned int lowest_cost_node_this_iteration = (-1);
      unsigned int lowest_cost_src_node_this_iteration = (-1);

      // iterate through all the solution points
      for(int j=0; j<solution_points_found; j++)
      {

         // find the node info for the node being examined (this lookup cannot fail)
         node_type::iterator itNode=graphNodes.find(mst_for_graph[j][DST_NODENUM_IDX]);
        
         // look through the edges of all the nodes in the solution so far
         for(edge_type::iterator itEdge=itNode->second->m_edges.begin(); itEdge != itNode->second->m_edges.end(); ++itEdge)
         {

            // don't look at connected nodes that are already in the solution (aka visited) OR
            //  nodes that aren't the right color for this tree
            if(isNodeVisited(itEdge->first) || getNodeColor(itEdge->first) != color)
            {
                   continue;
            }

             // now check if its the lowest cost so far in this iteration
            // if so, then record it
            if((itEdge->second) < lowest_cost_edge_this_iteration)
            {
               lowest_cost_edge_this_iteration = itEdge->second;
               lowest_cost_node_this_iteration = itEdge->first;
               lowest_cost_src_node_this_iteration = itNode->first;
            }
         }
      }

      if(lowest_cost_node_this_iteration == (-1))
      {
         break;
      }

      // add a newly found lowest node to the solution
      mst_for_graph[solution_points_found][DST_NODENUM_IDX] = lowest_cost_node_this_iteration;
      mst_for_graph[solution_points_found][SRC_NODENUM_IDX] = lowest_cost_src_node_this_iteration;
      mst_for_graph[solution_points_found][EDGEWEIGHT_IDX] = lowest_cost_edge_this_iteration;
      setNodeVisited(lowest_cost_node_this_iteration);

      // if the caller wanted the results, the record them
      if(pSolutionVector) pSolutionVector->push_back(mst_for_graph[solution_points_found][DST_NODENUM_IDX]);

      solution_points_found++;

   }

   if(printSolution)
   {

      // the MST solution is done, now print it out
      //
      std::cout << "----------- MST path ---------------" << std::endl;

      unsigned int mst_total_cost = 0;

      // iterate through the solution vector...
      for(i=1; i<solution_points_found; i++)
      {
         std::cout << "Node: " << mst_for_graph[i][SRC_NODENUM_IDX]  << "-->" <<  mst_for_graph[i][DST_NODENUM_IDX]
                   << " \tCost: " << mst_for_graph[i][EDGEWEIGHT_IDX] <<std::endl;
         mst_total_cost += mst_for_graph[i][EDGEWEIGHT_IDX];
      }

      // don't forget the total MST cost!
      std::cout << "\nTotal MST cost: " << mst_total_cost << std::endl;
   }

   return true;

}

void Graph::doDijkstra( unsigned int originNode, unsigned int destNode, std::list<unsigned int> *pathList, int &pathCost)
{
   bool validRouteFoundToDestination = false;
   unsigned int closedNodeNum = originNode;

   // initialize the outcome
   pathCost = 0;
   pathList->clear();

   // special case for origin == destination, just return
   if(originNode == destNode)
   {
      return;
   }

   // before starting, clean the nodes of computed values in case we're re-running the algorythm
   for(std::map<int, graphPoint* >::iterator itGraphNode = graphNodes.begin(); itGraphNode != graphNodes.end(); ++itGraphNode)
   {
      itGraphNode->second->cleanNode();
   }

   //
   // run the shortest path algorithm on the graph passed in
   //
   makeOriginNode(originNode);

   std::vector<unsigned int> openSet;    // a set of all "not yet visited" Nodes


   while(true)
   {
      int numNodesAddedToOpenSet = 0;
         
      //
      // add the connected nodes from this node to the open set (Step N+1)
      //
      std::map<int, graphPoint* >::iterator itGraphNode = graphNodes.find(closedNodeNum);

      unsigned int closedNodeCost = itGraphNode->second->m_totalCost;


      // get all the nodes connected to this one and 
      // adjust the costs and (from node) values of each connected node
      // itGraphEdge->first is the connected node number and itGraphEdge->second is the edge cost
      for(std::map<unsigned int, unsigned int>::iterator itGraphEdge= itGraphNode->second->m_edges.begin(); 
          itGraphEdge != itGraphNode->second->m_edges.end(); 
          ++itGraphEdge)
      {

         // first, mark this node as visited
         itGraphNode->second->setVisited();

         std::map<int, graphPoint* >::iterator itNextEdgeNode = graphNodes.find(itGraphEdge->first);

         if(itNextEdgeNode == graphNodes.end()) continue;  // no node actually exists 

         int i;
         for(i=0; i<openSet.size(); i++)
         {
            if(openSet[i] == itNextEdgeNode->second->m_nodeNumber) break;
         }

          // add it to the open set if it's not there already and it hasn't already been visited
         if(i == openSet.size() && (itNextEdgeNode->second->getVisited() == false))
         {
            openSet.push_back(itNextEdgeNode->second->m_nodeNumber);
            numNodesAddedToOpenSet++;
          }

          // now see if this is a lower cost path to this connected node
          if(itNextEdgeNode->second->getPointCost() == (-1) || 
             (closedNodeCost + itGraphEdge->second) <  itNextEdgeNode->second->getPointCost())
          {

              // new lower cost found
             itNextEdgeNode->second->setPointCost(static_cast<int>(closedNodeCost + itGraphEdge->second));
             itNextEdgeNode->second->m_viaNode = closedNodeNum;


          } 
       }

       //if we didn't add any new members to the openset and it's empty, we are done
       if ( (openSet.size() == 0) &&  (numNodesAddedToOpenSet == 0) )
       {
          break;
       }

       // Now, find the lowest cost member of the open set and make it the new closed set member

       unsigned int lowest_cost_node;
       unsigned int lowest_cost_index;

       for (int i=0; i<openSet.size(); i++)
       {
          int lowest_cost_val;

          std::map<int, graphPoint* >::iterator itGraphNode = graphNodes.find(openSet[i]);

          if((i==0) || (itGraphNode->second->m_totalCost < lowest_cost_val ))
          {
             lowest_cost_node = itGraphNode->second->m_nodeNumber;
             lowest_cost_val = itGraphNode->second->m_totalCost;
             lowest_cost_index = i;
          }
       }


       // this is the node that we'll be evaluating in the next iteration
       closedNodeNum = lowest_cost_node;

       // remove that node from the open set
       for(std::vector<unsigned int>::iterator vit = openSet.begin(); vit != openSet.end(); ++vit)
       {

          if (*vit == lowest_cost_node)
          {
             openSet.erase(vit);
             break;
          }
       }

        //if that node is the dest node, we have succeeded and we are done
       if (closedNodeNum == destNode)
       {
          validRouteFoundToDestination = true;
          break;
       }

       //else keep on truckin'
    }

   //
   // Now tell the user whether or not we've been able to find a route
   //   If so, print out the route
   //   else, well..there's not much we can do except break the bad news
   //
    if(validRouteFoundToDestination == true)
    {
       std::vector<unsigned int> routeSet;

       // now print out the route

       unsigned int routeNode = destNode;
       unsigned int lastRouteNode = routeNode;

       // reuse the openSet vector to record the final route
       pathCost = 0;

       while(true)
       {
          pathList->push_front(routeNode);

          // if we've added the orig to the route record, then we're done
          if(routeNode == originNode) break;

          // remember the "from" node for the cost computation
          lastRouteNode = routeNode;

          // find the "next" node
          routeNode = graphNodes.find(routeNode)->second->m_viaNode;

          // record this hop in the total cost of the path
          pathCost += getEdgeValue(routeNode, lastRouteNode);

          
       }
    }
    else
    {
       std::cout << "Could not find a route from " << originNode << " to " << destNode << std::endl;
       pathList->clear();  // no elements in the list
       pathCost = -1;
    }

}

inline std::ostream &right_shift_row(std::ostream &out, unsigned int row, bool printRow = false)
{
   
   if(printRow)
   {
      out << row;
   }

   out << "   ";

   // account for printing the row number
   if(printRow)
      if(row > 9) out<<"\b\b";
      else out<<"\b";
      
   // now print the spaces needed for each row for correct allignment
   for(int spaces=0; spaces<(row*2); spaces++)
   {
      out << " ";
   }
      return out;
}


ShortestPathAlgo::ShortestPathAlgo() : pathCost(-1)
{
   pathList = new std::list<unsigned int>;
} 

ShortestPathAlgo::~ShortestPathAlgo()
{
   delete pathList;
} 
 
// returns a count of the nodes
unsigned int ShortestPathAlgo::verticies(Graph &G)
{
   return G.getNodeCount();
}

// returns the cost of the path (or -1 if no path exists)
int ShortestPathAlgo::path_size( Graph &G, unsigned int originNode, unsigned int destNode )
{
   G.doDijkstra(originNode, destNode, pathList, pathCost);
   return pathCost;
}

// returns a list with the path
std::list<unsigned int> *ShortestPathAlgo::path( Graph &G, unsigned int originNode, unsigned int destNode)
{
   G.doDijkstra(originNode, destNode, pathList, pathCost);
   return pathList;
}

std::ostream &operator<< (std::ostream &cout, std::list<unsigned int> *path)
{
   unsigned routeLen = path->size();

   if(path->size())
   {
      cout << "==== Shortest Route has " << routeLen << " nodes" << std::endl;

      cout << "===== Shortest Route is : ";

      // There's a valid route, so iterate through the list and print each member...
      for(std::list<unsigned int>::iterator listIt = path->begin(); listIt != path->end(); ++listIt)
      {
         cout << *listIt << " ";
      }
   }
   else
   {
      cout << "No route found" << std::endl;
   }
      
   return cout;
}

// ========================================
// The HexGame class implemenatation
// ========================================
//
//

void hexGame::construct_board()
{
 
   m_pBoard = new Graph();
   int neighbor_array[3][2] = { { 0, 1}, {-1, 1}, {-1, 0}};
   
   // first create the basic hex board with no connections
   for(int row=0; row < m_boardsize; row++)
   {
      for(int col=0; col < m_boardsize; col++)
      {
         m_pBoard->addNode(getTileNumber(row,col));
      }
   }

   // now add all the edges
   // first create the basic graph node framework
   for(int row=0; row < m_boardsize; row++)
   {
      for(int nodeInRow=0; nodeInRow < m_boardsize; nodeInRow++)
      {
         for (int rowNeighbor = -1; rowNeighbor < 2; rowNeighbor++ )  //row above, this row, and next row
         {
            for (int neighbor_index = 0; neighbor_index < 2; neighbor_index++)  // neighbor to the left and right
            {

               // ignore the exception cases
               if(((rowNeighbor + row) < 0) || (rowNeighbor + row) > (m_boardsize - 1)) continue;
               if((( neighbor_array[rowNeighbor+1][neighbor_index] + nodeInRow ) < 0) || 
                  (( neighbor_array[rowNeighbor+1][neighbor_index] + nodeInRow ) > (getGameSize() - 1))) continue;
  
               // add the node to the graph
               m_pBoard->addEdge(((row * getGameSize()) + nodeInRow), 
                  (((rowNeighbor + row) * getGameSize()) + ((neighbor_array[rowNeighbor+1][neighbor_index] + nodeInRow ))),1);
            }
         }
      }
   }
}

//
// hexGame constructors
//
hexGame::hexGame(unsigned int boardsize):m_boardsize(boardsize)
{
   construct_board();
}

hexGame::hexGame()
{


   std::cout << "\n\nGame of Hex..." << std::endl;
   std::cout << "Choose a board size ";
   std::cin >> m_boardsize;

   // construct the game board
   construct_board();

}

hexGame::~hexGame()
{
   std::cout << "hexGame destructor called" << std::endl;
   delete m_pBoard;   // clean up class instance
}

inline unsigned int hexGame::getGameSize()
{
   return m_boardsize;
}

// zero based rows, col, and result
inline int hexGame::getTileNumber(unsigned int row, unsigned int col)
{
   if((row > m_boardsize-1) || ( col > m_boardsize-1))
   {
      return -1;
   }
   else
   {
      return ((row * m_boardsize) + col);
   }
}

// zero based
int hexGame::getTileRow(unsigned int nodeNumber)
{
   if(nodeNumber > (( m_boardsize *  m_boardsize)-1)) return -1;
   else return (nodeNumber / m_boardsize);
}

int hexGame::getTileColumn(unsigned int nodeNumber)
{
   if(nodeNumber > (( m_boardsize *  m_boardsize)-1)) return -1;
   else return (nodeNumber % m_boardsize );

}


nodecolor hexGame::getTileColor(unsigned int nodeNumber)
{
   return m_pBoard->getNodeColor(nodeNumber);
}


bool hexGame::isMoveLegal(unsigned int row, unsigned int col)
{
   if(m_pBoard->getNodeColor(getTileNumber(row, col)) == nodecolor::NONE) return true;
   else return false;
}

bool hexGame::makeGameMove(HexPlayer &p, unsigned int row, unsigned int col)
{
   bool ret = false;

   if(isMoveLegal(row,col))
   {
      m_pBoard->setNodeColor(getTileNumber(row, col), p.getPlayerColor());
      p.setLastMove(getTileNumber(row, col));
      ret = true;
   }

   return ret;

}


bool hexGame::makeGameMove(HexPlayer &p, unsigned int node)
{
   return makeGameMove(p, getTileRow(node) , getTileColumn(node));
}

void hexGame::forceTileColor(unsigned int node, playercolor color )
{
   m_pBoard->setNodeColor(node, color);
}


bool hexGame::checkForWinner(HexPlayer *player)
{

   // do we have a winner?
   // run the spanning tree on the user's last choice
   // if the RED player has MST nodes on north and south edges of the graph then they win
   // if the BLUE player has MST nodes on the east and west edges of the graph, then they win
   // else...PLAY ON

   std::vector< int > mst_solution_vec;
   bool ret = false;

   m_pBoard->doPrim( player->getLastMove(), 
      player->getPlayerColor(),
      &mst_solution_vec, false);

   // now see if we have a winner
   if(player->getPlayerColor() == nodecolor::BLUE)
   {
      // check if we have a solution that has both north and south edges
      if(mstIncludesRow(0, &mst_solution_vec) && 
         mstIncludesRow((getGameSize()-1), &mst_solution_vec))
      {
         ret = true;
      }
   }
   else
   {
      // check if we have a solution that has both east and west edges
      if(mstIncludesColumn(0, &mst_solution_vec) && 
         mstIncludesColumn((getGameSize()-1), &mst_solution_vec))
      {
         ret = true;
      }
   }

   return ret;

}

bool hexGame::mstIncludesColumn(unsigned int column, std::vector< int > *pMstVector)
{
   for(int i=0; i< pMstVector->size(); i++)
   {
      if(getTileColumn((*pMstVector)[i]) == column) return true;
   }

   return false;
}

bool hexGame::mstIncludesRow(unsigned int row, std::vector< int > *pMstVector)
{
   for(int i=0; i< pMstVector->size(); i++)
   {
      if(getTileRow((*pMstVector)[i]) == row) return true;
   }

   return false;
}

// on entry, the vector should be empty..but just in case...clear it out
bool hexGame::getAllFreeTiles(std::vector< unsigned int > &FreeTiles)
{
   FreeTiles.clear();

   for(int node=0; node < (getGameSize() * getGameSize()); node++)
   {
      if(getTileColor(node) == nodecolor::NONE)
      {
         FreeTiles.push_back(node);
      }
   }

   return true;
}




// overloaded output IO operator for a hex game board
std::ostream& operator<<(std::ostream &out, hexGame &game) 
{
   // print the graph
   // for hex graphs, this needs to be done row by row and indented

   out << "\n ============ Game of HEX (" << game.getGameSize() << " x " << game.getGameSize() << ")===============" << std::endl;
   out << " \e[1;31mRED plays side-side\e[0m, and \e[1;36mBLUE plays up-down\e[0m\n\n";


   // hex graphs have a dimension greater than 0
   if(game.getGameSize() > 0)
   {
      // this is a hex graph
      right_shift_row(out, 1);
      for(int header=0; header<game.getGameSize(); header++) 
      {
         if(!(header<9)) out <<"\b";  // adjust for the size of the number
         out << (header+1) << "   ";

      }
      out << "\n\n";

      // must do two iterations per row
      for(int row=0; row<game.getGameSize(); row++)
      {

         // we need to do two passes for each row, one for the nodes, and one for the interconnects
         for(int rowIter=0; rowIter<2; rowIter++)
         {
            // print the row number on the node row and allign correctly
            right_shift_row(out, row+1, (rowIter ? false : true));

            // do once for each column in the row
            for(int col=0; col<game.getGameSize();col++)
            {
               unsigned int nodeNumber = game.getTileNumber(row, col);

               if(game.m_pBoard->isNodeInGraph(nodeNumber))
               {
                  if(!rowIter)
                  {
                     // we're printing the even row 
                     // first the node
                     out << ((game.getTileColor(nodeNumber) == (nodecolor::NONE)) ? 
                        "*" : (game.getTileColor(nodeNumber) == (nodecolor::RED)) ? "\e[1;31mX\e[0m" : "\e[1;36mO\e[0m");
                     
                     // then the neighor on this row
                     int neighborNodeNumber = game.getTileNumber(row, col+1);

                     // out << "node " << nodeNumber << "neighbor" << neighborNodeNumber << std::endl;
                     if(neighborNodeNumber != -1)
                     {
                        // the UGLIEST nested tri-graph in the history of ugly trigraphs!!!
                        out << ((game.m_pBoard->hasEdge(nodeNumber, neighborNodeNumber)) ? 
                           (game.getTileColor(nodeNumber) == game.getTileColor(neighborNodeNumber) ? 
                              (game.getTileColor(nodeNumber) == nodecolor::NONE ? 
                                 "\e[1;33m - \e[0m" : 
                                 (game.getTileColor(nodeNumber) == nodecolor::RED ? 
                                    "\e[1;31m - \e[0m" : "\e[1;36m - \e[0m")) : "\e[1;33m - \e[0m" ) : " " );
                     }
                                                                   
                  }
                  else
                  {
                     // print the between rows
                     int neighborNodeNumber = game.getTileNumber(row+1, col-1);
                     if(neighborNodeNumber != -1)
                     {
                        out << ((game.m_pBoard->hasEdge(nodeNumber, neighborNodeNumber)) ? 
                           (game.getTileColor(nodeNumber) == game.getTileColor(neighborNodeNumber) ? 
                              (game.getTileColor(nodeNumber) == nodecolor::NONE ? 
                                 "\e[1;33m/ \e[0m" : 
                                 (game.getTileColor(nodeNumber) == nodecolor::RED ? 
                                    "\e[1;31m/ \e[0m" : "\e[1;36m/ \e[0m")) : "\e[1;33m/ \e[0m" ) : " " );
                     }
                     else
                     {
                        out << " ";
                     }

                     // print the between rows
                     neighborNodeNumber = game.getTileNumber(row+1, col);
                     if(neighborNodeNumber != -1)
                     {
                        out << ((game.m_pBoard->hasEdge(nodeNumber, neighborNodeNumber)) ? 
                           (game.getTileColor(nodeNumber) == game.getTileColor(neighborNodeNumber) ? 
                              (game.getTileColor(nodeNumber) == nodecolor::NONE ? 
                                 "\e[1;33m\\ \e[0m" : 
                                 (game.getTileColor(nodeNumber) == nodecolor::RED ? 
                                    "\e[1;31m\\ \e[0m" : "\e[1;36m\\ \e[0m")) : "\e[1;33m\\ \e[0m" ) : " " );
                      }


                  }
               }
            }
            out << "\n";
         }
      }
   }

   return out;
}




// ===================================
// the Player class implementation
// ===================================
//
Player::Player(){}

Player::Player(std::string name)
{
   m_name = name;
}

void Player::setPlayerIndex(unsigned int player_index)
{
   index = player_index;
}

unsigned int Player::getPlayerIndex()
{
   return index;
}

// ===================================
// the HexPlayer class implementation
// ===================================
//

HexPlayer::HexPlayer(hexGame &game, std::string name):m_game(game), Player(name){}

playercolor HexPlayer::getPlayerColor()
{
   return m_color;
}

void HexPlayer::setPlayerColor(playercolor color)
{
   std::cout << m_name << " player is " << color << std::endl;
   m_color = color;
}

playertype HexPlayer::getPlayerType()
{
   return m_type;
}

void HexPlayer::setPlayerType(playertype type)
{
   m_type = type;
}

void HexPlayer::setLastMove(unsigned int nodenumber)
{
   m_lastMove = nodenumber;
}

unsigned int HexPlayer::getLastMove()
{
   return m_lastMove;
}


// ========================================
// The HumanHexPlayer class implemenatation
// ========================================
//
HumanHexPlayer::HumanHexPlayer(hexGame &game):HexPlayer(game, "human"){}

void HumanHexPlayer::makeMove()                             // makes a move appropriate for a human player
{
   unsigned int choice_row, choice_column;

   do
   {
      std::cout << m_color << ", choose a column: ";
      std::cin >> choice_column;

      if((choice_column == 0) || (choice_column > m_game.getGameSize()))
      {
         std::cout << "you must choose a column between 1 and " << (m_game.getGameSize()) << std::endl;
         continue;
      }

      break;  // we're done with the column choice when we get a good response from the user

   }while(true);

   do
   {
      std::cout << m_color << ", choose a row: ";
      std::cin >> choice_row;

      if((choice_row == 0) || (choice_row > m_game.getGameSize()))
      {
         std::cout << "you must choose a row between 1 and " << (m_game.getGameSize()) << std::endl;
         continue;
      }

      // check to see if the node already has a color
      if(m_game.getTileColor(m_game.getTileNumber(choice_row-1, choice_column-1)) != nodecolor::NONE)
      {
         std::cout << "Your tile choice has already been taken, please try again" << std::endl;
         continue;
                 
      }

      // either try again, or we're done
      break;

   }while(true);

   // make the game move for this player
   m_game.makeGameMove(*this, (choice_row-1), (choice_column-1));


}


// ===========================================
// The ComputerHexPlayer class implemenatation
// ===========================================
ComputerHexPlayer::ComputerHexPlayer(hexGame &game):HexPlayer(game, "computer")
{
   mc_sim = new monteCarlo(m_game);
}

void ComputerHexPlayer::makeMove()                            // makes a move appropriate for a computer player
{
   unsigned int next_move;

   // call the monte-carlo simulation for this player
   next_move = (*mc_sim)(*this);

   // and make the move that the computer has chosen
   m_game.makeGameMove(*this, m_game.getTileRow(next_move), m_game.getTileColumn(next_move));


}


//
// main entry point to the code
//
int main()
{

    unsigned int board_size;

    /* initialize random seed: */
    srand (time(NULL));

    const char print_graph_entry = 'n';

    // create a new game
    hexGame hexgame;

    // our two players
    std::array< HexPlayer *, 2 > hex_players;

    // First, let player one choose the color he wants
    HexPlayer *p1 = new HumanHexPlayer(hexgame);

    char players_color_choice;

    std::cout << "Choose a color \e[1;31m(R) for Red\e[0m, \e[1;36m(B) for Blue\e[0m: ";
    do{
       std::cin >> players_color_choice;
       if((players_color_choice == 'R') || (players_color_choice == 'r'))
       {
          p1->setPlayerColor(playercolor::RED);
          p1->setPlayerIndex(1);
          break;
       }
       else if ((players_color_choice == 'B') || (players_color_choice == 'b'))
       {
          p1->setPlayerColor(playercolor::BLUE);
          p1->setPlayerIndex(0);             
          break;
       }
       else
       {
          std::cout << "Please choose either \e[1;31m(R) for Red\e[0m or (\e[1;36m(B) for Blue\e[0m..." << std::endl;
       }
    }while(true);

    // store the player instance
    hex_players[p1->getPlayerIndex()] = p1;

          
    // Second player gets the color that's left...
    HexPlayer *p2 = new ComputerHexPlayer(hexgame);

    p2->setPlayerIndex( p1->getPlayerIndex() ? 0 : 1);
    p2->setPlayerColor(p1->getPlayerColor() == playercolor::RED ? playercolor::BLUE : playercolor::RED);
    hex_players[p2->getPlayerIndex()] = p2;

    
    // print out the board
    std::cout << hexgame << std::endl;

    // keep playing until there's a winner
    while(true)
    {

       for(int index=0; index<2; index++)
       {

          // allow each player to make their move...blue goes first
          hex_players[index]->makeMove();

          if(hexgame.checkForWinner(hex_players[index]))
          {
             std::cout << hex_players[index]->getPlayerColor() << " Wins!!" << std::endl;
             exit(0);
          }
          
          // print out the board
          std::cout << hexgame << std::endl;
          
       }
    }

}

