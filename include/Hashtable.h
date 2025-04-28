#ifndef HASHTABLE_H
#define HASHTABLE_H

#include <cstdint>
#include <vector>
#include <iostream>
#include <glm/vec2.hpp>
#include "sharedVariables.h"

class HashTable {
	
public:

	HashTable(unsigned int numberOfParticles);

	std::vector<unsigned int> computeNeighboringCellsKeys(glm::vec2& position);
	std::vector<unsigned int>& getCell(unsigned int key) { return particleTable[key]; };

	/*
	
		Compute neighbors of all particles
	*/

	/*
	void computeAllNeighbors(std::vector<glm::vec2>& particlePositions, std::vector<std::vector<unsigned int>>& neighbors);
	
	*/
	void insertInCell(glm::vec2& p, unsigned int index);

	/*
		initializeTable the initial form of the table is a n x n square pool of particles
		....
		....
		....
		....
	*/

	void initializeTable(std::vector<glm::vec2>& positions);

	void clearTable();

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
	std::vector<std::vector<unsigned int>> particleTable;
	const int P1 = 73856093;
	const int P2 = 19349663;
	int numParticles;

	// The spacing between cell
	float cellSize {0.2f};

	struct 
	{
		unsigned int horizontalCellsCount{40};
		unsigned int verticalCellsCount{40};
	} boundedTableInfo;

	

	/*
		cellDisplacements represent the surrounding cells plus the current cell
	*/
	int cellDisplacements[9][2] = { {0, 0}, {-1, 0}, {1, 0}, {0, -1}, {0, 1}, {-1, 1}, {1, 1}, {-1, -1}, {1, -1} };


	int tableHash(glm::vec2& cellPos);

	glm::vec2 getCellPosition(glm::vec2& position);


	int computeKey(glm::vec2& position);
	/*
		Compute the neighbors of a single particle
	*/
	void getNeighbors(std::vector<glm::vec2>& positions, int currentIndex, std::vector<std::vector<unsigned int>>& nbs);

};


#endif // !HASHTABLE_H
