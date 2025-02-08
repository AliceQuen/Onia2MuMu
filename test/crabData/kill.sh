#!/bin/bash

# 检查是否提供了参数
if [ -z "$1" ]; then
    echo "Usage: $0 <task_number>"
    exit 1
fi

# 获取所有匹配的crab任务的列表
tasks=$(ls -d crab_crab3_"$1"_*)

# 遍历所有任务并kill
for task in $tasks; do
    if [ -d "$task" ]; then
        echo "Killing task: $task"
        crab kill -d "$task"
    else
        echo "No valid CRAB project directories found for crab3_$1."
    fi
done

echo "All crab3_$1 tasks have been killed."