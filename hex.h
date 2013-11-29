#ifndef HEX_H
#define HEX_H


// forward class declarations
class graphPoint;
#include "graph.h"
class Player;

// the playertype definition
enum class playertype { 
   HUMAN, COMPUTER };

typedef nodecolor playercolor;    // the node colors are also the player colors
typedef unsigned int movevalue;   // the resultant value of a generic move (range is 0-worst to 100-best)


class minimax {
 private:
   Graph tree;                    // the minimax tree instance to use in this class.

 public:
   minimax();                                            // class constructor
   void addNode(unsigned int NodeNumber);
   void connectToNode(unsigned int srcNodeNumber, unsigned int destNodeNumber);  // connects two nodes together
   unsigned int operator()(unsigned int topNodeNumber);  // run the algorigthm itself from the top node, return the value found
};


class hexGame {

 private:
   const unsigned int boardsize;  // the size of the hex playing board
   const Graph *pGameboard;       // a pointer to the playing board itself

 public:
   hexGame();
   ~hexGame();    // destructor to free up allocated resources 
   bool makeMove(playercolor color, unsigned int col, unsigned int row);    // returns true if move was good
   bool isMoveLegal(unsigned int col, unsigned int row);  // returns true if move is legal
   unsigned int evaluateMove(unsigned int col, unsigned int row, Player &p); // returns value of move for player p

};


// A generic player.  Type of game unknown.
class Player 
{
 protected:
   const unsigned int index;                        // the assigned index for the player (0-none)
   playertype type;                                 // the type of player (human, computer)

 public:
   virtual void makeMove() = 0;                // pure virtual function for making a move
   unsigned int getPlayerIndex();
   void setPlayerIndex(unsigned int player_index);
   
   void setPlayerType(playertype type);
   playertype getPlayerType();
};

// A hex game player
class HexPlayer:public Player 
{
protected:
   const playercolor playerColor;                  // playercolor::RED or playercolor::BLUE
   const hexGame &GameBoard;                       // a reference to the gameboard being used in this game

public:
   HexPlayer(const playercolor color, const hexGame &board );
   playercolor getPlayerColor();
   playercolor setPlayerColor(playercolor color);
   virtual void makeMove() = 0;                    // pure virtual move function
   
};

// A computer hex game player
class ComputerHexPlayer:public HexPlayer
{
 private:
   HexGame &gameboard;

 public:
   ComputerHexPlayer(const hexGame &board );
   virtual void makeMove();                        // makes a move appropriate for a computer player
   
};

// a human hex game player
class HumanHexPlayer:public HexPlayer
{

 public:
   HumanHexPlayer(const hexGame &board);
   virtual void makeMove();                        // makes a move appropriate for a human player
   
};

#endif
