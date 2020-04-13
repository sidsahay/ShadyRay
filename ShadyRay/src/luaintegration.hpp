//
// Created by walksbynight on 12/4/18.
//

#ifndef SHADYRAY_LUAINTEGRATION_HPP
#define SHADYRAY_LUAINTEGRATION_HPP

#include <iostream>
#include <vector>
#include <cstring>

extern "C" {
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
}

typedef lua_State* LUA_state;


template <typename T>
T LUA_return(LUA_state L) {
    return 0;
}

template <>
const char * LUA_return<const char *>(LUA_state L) {
    return lua_tostring(L, -1);
}

template <>
int LUA_return<int>(LUA_state L) {
    return (int)lua_tointeger(L, -1);
}

template <>
double LUA_return<double>(LUA_state L) {
    return lua_tonumber(L, -1);
}

template <typename T>
T LUA_get_table_field(LUA_state L, const char *table_accessor) {
    char str[128] = {0};
    strcpy(str, table_accessor);

    std::vector<char*> tokens;
    char *token = strtok(str, ".");
    while (token != nullptr) {
        tokens.push_back(token);
        token = strtok(nullptr, ".");
    }

    for (int i = 0; i < tokens.size(); i++) {
        if (i == 0) {
            lua_getglobal(L, tokens[0]);
        }
        else {
            lua_pushstring(L, tokens[i]);
            lua_gettable(L, -2);
        }
    }

    auto ans = LUA_return<T>(L);
    lua_pop(L, tokens.size());
    return ans;
}

template <typename T>
T LUA_get_table_field_local(LUA_state L, const char *table_accessor) {
    char str[128] = {0};
    strcpy(str, table_accessor);

    std::vector<char*> tokens;
    char *token = strtok(str, ".");
    while (token != nullptr) {
        tokens.push_back(token);
        token = strtok(nullptr, ".");
    }

    for (int i = 0; i < tokens.size(); i++) {
        if (i == 0) {
            //lua_getglobal(L, tokens[0]);
        }
        else {
            lua_pushstring(L, tokens[i]);
            lua_gettable(L, -2);
        }
    }

    auto ans = LUA_return<T>(L);
    lua_pop(L, tokens.size() - 1);
    return ans;
}

LUA_state LUA_initialize(const char *file_name) {
    LUA_state L = luaL_newstate();
    luaL_openlibs(L);

    if (luaL_loadfile(L, file_name) != LUA_OK) {
        std::cout << "Lua: " << lua_tostring(L, -1) << std::endl;
    }

    return L;
}

void LUA_execute(LUA_state L) {
    if (lua_pcall(L, 0, 0, 0) != LUA_OK) {
        std::cout << "Lua Exec: " << lua_tostring(L, -1) << std::endl;
    }
}

void LUA_register_function(LUA_state L, const char *name, int(*foo)(LUA_state)) {
    lua_register(L, name, foo);
}

void LUA_close(LUA_state L) {
    lua_close(L);
}

class LuaContext {
public:
    explicit LuaContext(LUA_state L);
    explicit LuaContext(const char *file_name);

    template <typename T> T GetValue();
    template <typename T> T GetTableField(const char *table_accessor);
    template <typename T> T GetTableFieldLocal(const char *table_accessor);

    void RegisterFunction(const char *name, int(*foo)(LUA_state L));

    void Execute(int num_args, int num_return);

    void Pop(int n);

    LUA_state GetState();

    ~LuaContext();

private:
    bool is_external_context = false;
    LUA_state L;
};

LuaContext::LuaContext(LUA_state L) : L(L), is_external_context(true){
}

LuaContext::LuaContext(const char *file_name) : is_external_context(false) {
    L = luaL_newstate();
    luaL_openlibs(L);

    if (luaL_loadfile(L, file_name) != LUA_OK) {
        std::cout << "Lua: " << lua_tostring(L, -1) << std::endl;
    }
}

template<typename T>
T LuaContext::GetValue() {
    return nullptr;
}

template <>
int LuaContext::GetValue<int>() {
    return (int)lua_tointeger(L, -1);
}

template <>
double LuaContext::GetValue<double>() {
    return lua_tonumber(L, -1);
}

template <>
const char* LuaContext::GetValue<const char*>() {
    return lua_tostring(L, -1);
}

template<typename T>
T LuaContext::GetTableField(const char *table_accessor) {
    char str[128] = {0};
    strcpy(str, table_accessor);

    std::vector<char*> tokens;
    char *token = strtok(str, ".");
    while (token != nullptr) {
        tokens.push_back(token);
        token = strtok(nullptr, ".");
    }

    for (int i = 0; i < tokens.size(); i++) {
        if (i == 0) {
            lua_getglobal(L, tokens[0]);
        }
        else {
            lua_pushstring(L, tokens[i]);
            lua_gettable(L, -2);
        }
    }

    auto ans = LUA_return<T>(L);
    lua_pop(L, tokens.size());
    return ans;
}

template<typename T>
T LuaContext::GetTableFieldLocal(const char *table_accessor) {
    char str[128] = {0};
    strcpy(str, table_accessor);

    std::vector<char*> tokens;
    char *token = strtok(str, ".");
    while (token != nullptr) {
        tokens.push_back(token);
        token = strtok(nullptr, ".");
    }

    for (int i = 0; i < tokens.size(); ++i) {
        if (i == 0) {
            //lua_getglobal(L, tokens[0]);
        }
        else {
            lua_pushstring(L, tokens[i]);
            lua_gettable(L, -2);
        }
    }

    auto ans = LUA_return<T>(L);
    lua_pop(L, tokens.size() - 1);
    return ans;
}

void LuaContext::RegisterFunction(const char *name, int (*foo)(LUA_state)) {
    LUA_register_function(L, name, foo);
}

void LuaContext::Execute(int num_args, int num_return) {
    if (lua_pcall(L, num_args,  num_return, 0) != LUA_OK) {
        std::cout << "Lua Exec: " << lua_tostring(L, -1) << std::endl;
    }
}

LUA_state LuaContext::GetState() {
    return L;
}

LuaContext::~LuaContext() {
    if (!is_external_context) {
        LUA_close(L);
    }
}

void LuaContext::Pop(int n) {
    lua_pop(L, n);
}

#endif //SHADYRAY_LUAINTEGRATION_HPP
