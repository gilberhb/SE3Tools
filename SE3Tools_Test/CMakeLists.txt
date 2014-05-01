SET(SOURCES	SE3Tools_Test.cpp)
include_directories(../SE3Tools)
add_executable(SE3Tools_Test	${SOURCES})
target_link_libraries(SE3Tools_Test	SE3Tools)
add_custom_command(TARGET SE3Tools_Test	POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy_if_different
	"$<TARGET_FILE:SE3Tools>"
	"$<TARGET_FILE_DIR:SE3Tools_Test>")

