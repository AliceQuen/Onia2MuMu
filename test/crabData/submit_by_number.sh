#!/bin/bash

# 检查是否提供了参数
if [ -z "$1" ]; then
    echo "Usage: $0 <0-7>"
    exit 1
fi

# 检查参数是否在0到7之间
if [[ "$1" -lt 0 || "$1" -gt 7 ]]; then
    echo "Parameter must be between 0 and 7."
    exit 1
fi

# 根据参数选择提交的数据集
for pyfile in crab3_"$1"_*_MINIAOD.py
do
    # 检查文件是否存在
    if [[ -f "$pyfile" ]]; then
        # 提交文件
        rm -r "./crab_${pyfile%.py}"
        echo "Submitting $pyfile"
        crab --quiet submit "$pyfile"
    else
        echo "No crab3_$1_*_MINIAOD.py files found."
    fi
done