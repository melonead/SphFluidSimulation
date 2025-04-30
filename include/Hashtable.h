#ifndef HASHTABLE_H
#define HASHTABLE_H

#include <cstdint>
#include <vector>
#include <iostream>
#include <glm/vec2.hpp>
#include "sharedVariables.h"


/*
		Query for the neighbors of a point
	*/
struct NeighborQuery
{
	std::vector<unsigned int> neighborBucket;
	unsigned int size{0};
};


class HashTable {
	
public:
	HashTable(unsigned int numberOfParticles);

	
	/*
		Get the IDs of all the neighboring particles.
	*/

	NeighborQuery getNeighborIDs(glm::vec2& position);
	

	void createTable(std::vector<glm::vec2>& positions);
	
	float getCellSize() {return cellSize; };
	unsigned int getHorizontalCellsCount() {return boundedTableInfo.horizontalCellsCount; };
	unsigned int getVerticalCellCount() {return boundedTableInfo.verticalCellsCount; };
	
private:
	/*
	particleTable represents a grid for faster computation of neighbors of a particles,
	To get the neighbors of a particle in a cell we need to only search the surrounding
	8 cell plus the current cell.
	
	A particle is represented by 1 number: 
	(1). Index in the position vector/list

	*/
	const int P1 = 73856093;
	const int P2 = 19349663;
	
	int numParticles;

	struct 
	{
		unsigned int horizontalCellsCount{40};
		unsigned int verticalCellsCount{40};
	} boundedTableInfo;
	
	/*
	Particle IDs, these do not change throught out the simulation.
	*/
	std::vector<unsigned int> particleIDs;
	/*
	Contain the number of particles in each cell. 
	By the end of the create table function, grid is going to tell us where the first
	of the particles in a particular cell are in the sortedParticleIDs.
	*/
	std::vector<unsigned int> grid;

	/*
	These ID are sorted in the sense that those that are in the same cell are next to each other.
	They change from frame to frame.
		*/
	std::vector<unsigned int> sortedParticleIDs;

	/*
		Since by the end of createTable, grid has transformed from telling us the number of particles in a particular
		cell to telling us where the first is, we need to store the count. particleCout does this.
	*/
	std::vector<unsigned int> particleCounts;

	NeighborQuery query;

	


	/*
	cellDisplacements represent the surrounding cells plus the current cell
	*/
	int cellDisplacements[9][2] = { {0, 0}, {-1, 0}, {1, 0}, {0, -1}, {0, 1}, {-1, 1}, {1, 1}, {-1, -1}, {1, -1} };

	int tableHash(glm::vec2& cellPos);

	glm::vec2 getCellPosition(glm::vec2& position);


	int computeKey(glm::vec2& position);
};


#endif // !HASHTABLE_H
