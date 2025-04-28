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
	unsigned numCells = getHorizontalCellsCount() * getVerticalCellCount();

	// Add an extra cell so that we don't index outside the array while computing
	// partial sums.
	grid = std::vector<unsigned int>(numCells + 1, 0);
	particleIDs = std::vector<unsigned int>(numberOfParticles, 0);
	sortedParticleIDs = std::vector<unsigned int>(numberOfParticles, 0);
	particleCounts = std::vector<unsigned int>(numCells, 0);

	// create particles IDs
	for (unsigned int i = 0; i < numberOfParticles; i++)
	{
		particleIDs[i] = i;
	}
}



std::vector<unsigned int> HashTable::getNeighborIDs(glm::vec2& position)
{
	// REVISIT: Perhaps make this a member variable?
	std::vector<unsigned int> ids;
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
		// Go to table and find where to start looking
		unsigned int start = grid[key];
		// Go to particleCount and find how much to advance start at start
		unsigned int size = particleCounts[key];
		for (unsigned int i = 0; i < size; i++)
		{
			unsigned int id = sortedParticleIDs[start + i];
			ids.push_back(id);
		}
	}

	return ids;
}

void HashTable::createTable(std::vector<glm::vec2>& positions)
{
	unsigned int size = positions.size();
	// Compute cell particle count
	for(unsigned int i = 0; i < size; i++)
	{
		glm::vec2& position = positions[i];
		unsigned int key = computeKey(position);
		grid[key]++;
		// This apparent duplication is necessary because cellParticlesCout is goint to change.
		particleCounts[key]++;
	}
	// compute partial sums
	unsigned int cellsNumber = getVerticalCellCount() * getHorizontalCellsCount();
	for(unsigned int i = 0; i < cellsNumber; i++)
	{
		grid[i + 1]  += grid[i];
	}
	// Sort the particle ID, particles in the same cell will be next to each other
	
	for(unsigned int i = 0; i < size; i++)
	{
		glm::vec2& position = positions[i];
		unsigned int particleID = particleIDs[i];

		unsigned int key = computeKey(position);

		grid[key]--;
		sortedParticleIDs[grid[key]] = particleID;
	}

	int k = 0;

}

// hash: compute hash key from the position of a particle
int HashTable::tableHash(glm::vec2& cellPos)
{
	int x = cellPos.x;
	int y = cellPos.y;

	int n = (getHorizontalCellsCount() * getVerticalCellCount());

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

int HashTable::computeKey(glm::vec2& position)
{
	glm::vec2 cellPos = getCellPosition(position);
	int key = tableHash(cellPos);

	return key;
}

/*
void HashTable::computeAllNeighbors(std::vector<glm::vec2>& particlePositions, std::vector<std::vector<unsigned int>>& neighbors) {
	for (unsigned int i = 0; i < numParticles; i++) {
		getNeighbors(particlePositions, i, neighbors);
	}
}
*/