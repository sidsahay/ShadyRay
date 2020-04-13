//
// Created by walksbynight on 23/3/18.
//

#ifndef SHADYRAY_MEMORYPOOL_HPP
#define SHADYRAY_MEMORYPOOL_HPP

#include <vector>

#define POOL_ALLOC(pool_obj, type, ctor) (new ((void*)((pool_obj).GetNewPointer(sizeof(type)))) ctor)
#define POOL_CLEAR(pool_obj) ((pool_obj)->Clear())

#define WITH_POOL(pool_name, size) for (MemoryPool pool_name(size); pool_name.Once();)

class ScopeDefinition {
public:
    ScopeDefinition();

    bool Once();

private:
    bool is_once;
};

class MemoryPool : public ScopeDefinition {
public:
    explicit MemoryPool(unsigned int pool_size);

    unsigned char *GetNewPointer(unsigned int object_size);

    void Clear();

    ~MemoryPool();

protected:
    std::vector<unsigned char *> pool;
    unsigned int pool_size;
    unsigned int current_free_index;
    unsigned int current_free_pool;
};

#endif //SHADYRAY_MEMORYPOOL_HPP
