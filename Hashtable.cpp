#include "Hashtable.h"
#include <cstdint>
#include <vector>
#include <iostream>
#include <glm/glm.hpp>
#include <cassert>
#include "sharedVariables.h"

HashTable::HashTable(unsigned int numberOfParticles)
	: numParticles(numberOfParticles)

{
	unsigned int numCells = getHorizontalCellsCount() * getVerticalCellCount();
	particleTable = std::vector<std::vector<unsigned int>>(numCells);
	for (unsigned int i = 0; i < (unsigned int) numCells; i++) {
		particleTable[i].reserve(10);
	}
	
}


void HashTable::initializeTable(std::vector<glm::vec2>& positions) {
    
	for(unsigned int i = 0; i < numParticles; i++) {
		insertInCell(positions[i], i);
	}
}

// hash: compute hash key from the position of a particle
int HashTable::tableHash(glm::vec2& cellPos)
{
	int x = cellPos.x;
	int y = cellPos.y;

	int n = (int)(getHorizontalCellsCount() * getVerticalCellCount());
	
	return abs(((x * P1) + ( y * P2)) % n);
}

// getCellPosition: convert world position to cell position
glm::vec2 HashTable::getCellPosition (glm::vec2& position)
{
	// discretize particle position from float to int
	int x = floor(position.x / radiusOfInfluence);
	int y = floor(position.y / radiusOfInfluence);
	return glm::vec2(x, y);
}

// vector<int> (&gr)[10]
// insertInMap: inserts the particle index into the hashtable
void HashTable::insertInCell(glm::vec2& p, unsigned int index)
{
	glm::vec2 cellPos = getCellPosition(p);
	int key = tableHash(cellPos);
	// std::cout << "key: " << key << " position: " << "{ " << p.x << ", " << p.y << " }" << std::endl;
	particleTable[key].push_back(index);
}

int HashTable::computeKey(glm::vec2& position)
{
	glm::vec2 cellPos = getCellPosition(position);
	int key = tableHash(cellPos);
	return key;
}

void HashTable::clearTable()
{
	for (int i = 0; i < (getHorizontalCellsCount() * getVerticalCellCount()); i++)
	{
		particleTable[i].clear();
	}
}

/*
void HashTable::getNeighbors(std::vector<glm::vec2>& positions, int currentIndex, std::vector<std::vector<unsigned int>>& nbs)
{
	glm::vec2 cellPos;
	glm::vec2 currentCellPos = getCellPosition(positions[currentIndex]);
	int key;
	for (int i = 0; i < 9; i++)
	{
		cellPos.x = currentCellPos.x + cellDisplacements[i][0];
		cellPos.y = currentCellPos.y + cellDisplacements[i][1];
		key = tableHash(cellPos);
		int size = particleTable[key].size();
		if (size == 0)
			continue;
		for (int j = 0; j < size; j++)
		{
			unsigned int nbIndex = particleTable[key][j];
			float xDist = positions[nbIndex].x - positions[currentIndex].x;
			float yDist = positions[nbIndex].y - positions[currentIndex].y;
			if ((xDist * xDist + yDist * yDist) < (radiusOfInfluence * radiusOfInfluence))
				nbs[currentIndex].push_back(nbIndex);

		}

	}
}
*/


std::vector<unsigned int> HashTable::computeNeighboringCellsKeys(glm::vec2& position) {
	std::vector<unsigned int>  keys;
	unsigned int neighborCellsCount = 9;
	glm::vec2 currentCellPos = getCellPosition(position);
	glm::vec2 neighborCellPos;
	unsigned int key;
	for (unsigned int i = 0; i < neighborCellsCount; i++) {
		// REVISIT: Do the cell displacements need to be scaled by
		// the cell size?
		neighborCellPos.x = currentCellPos.x + cellDisplacements[i][0];
		neighborCellPos.y = currentCellPos.y + cellDisplacements[i][1];
		key = tableHash(neighborCellPos);
		keys.push_back(key);
	}
	return keys;
}

/*
void HashTable::computeAllNeighbors(std::vector<glm::vec2>& particlePositions, std::vector<std::vector<unsigned int>>& neighbors) {
	for (unsigned int i = 0; i < numParticles; i++) {
		getNeighbors(particlePositions, i, neighbors);
	}
}
*/