SET(SOURCES	SE3.cpp
		SO3.cpp
		stdafx.cpp)
SET(HEADERS	SE3.h
		SO3.h
		stdafx.h)

add_library(SE3Tools SHARED ${SOURCES} ${HEADERS})
set_target_properties(SE3Tools PROPERTIES DEFINE_SYMBOL "BUILD_SHARED_LIB")
set_target_properties(SE3Tools PROPERTIES DEBUG_POSTFIX "d")