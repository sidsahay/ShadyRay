//
// Created by walksbynight on 23/3/18.
//

#include "memorypool.hpp"

ScopeDefinition::ScopeDefinition() : is_once(true) {
}

bool ScopeDefinition::Once() {
    if (is_once) {
        is_once = false;
        return true;
    } else {
        return false;
    }
}

MemoryPool::MemoryPool(unsigned int pool_size) : pool_size(pool_size), current_free_index(0), current_free_pool(0) {
    pool.push_back(new unsigned char[pool_size]);
}

MemoryPool::~MemoryPool() {
    for (auto ptr : pool) {
        delete[] ptr;
    }

    pool.clear();
}

unsigned char *MemoryPool::GetNewPointer(unsigned int object_size) {
    if ((current_free_index + object_size) < (pool_size - 1)) {
        unsigned char *ptr = pool[current_free_pool] + current_free_index;
        current_free_index += object_size;
        return ptr;
    } else {
        pool.push_back(new unsigned char[pool_size]);
        current_free_pool++;
        unsigned char *ptr = pool[current_free_pool];
        current_free_index = 0;
        return ptr;
    }
}

void MemoryPool::Clear() {
    for (auto ptr : pool) {
        delete[] ptr;
    }

    pool.clear();

    pool.push_back(new unsigned char[pool_size]);
    current_free_index = 0;
    current_free_pool = 0;
}
