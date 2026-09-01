#!/usr/bin/env python3
import re
from pathlib import Path

def extract_branch_variables(cc_file):
    """从 beginJob() 中提取所有 Branch 定义的变量"""
    with open(cc_file, 'r', encoding='utf-8') as f:
        content = f.read()
    
    begin_job_pattern = r'void MultiLepPAT::beginJob\(\)[^}]+?\} // begin Job'
    begin_job_match = re.search(begin_job_pattern, content, re.DOTALL)
    
    if not begin_job_match:
        return set()
    
    begin_job_content = begin_job_match.group(0)
    branch_pattern = r'X_One_Tree_->Branch\("([^"]+)",\s*&(\w+)'
    matches = re.findall(branch_pattern, begin_job_content)
    
    return {var for name, var in matches}

def extract_reset_variables(cc_file):
    """从 resetEventVariables() 和 resetCandidateVariables() 中提取所有被重置的变量"""
    with open(cc_file, 'r', encoding='utf-8') as f:
        content = f.read()

    all_reset_vars = set()

    for func_name in ('resetEventVariables', 'resetCandidateVariables', 'resetVariables'):
        # 使用花括号匹配来处理嵌套的 if/else 块
        pattern = r'void MultiLepPAT::' + func_name + r'\s*\(\s*\)\s*\{'
        match = re.search(pattern, content)
        if not match:
            continue

        start = match.end()
        brace_count = 1
        end = start
        for i in range(start, len(content)):
            if content[i] == '{':
                brace_count += 1
            elif content[i] == '}':
                brace_count -= 1
                if brace_count == 0:
                    end = i
                    break

        reset_content = content[start:end]
        var_pattern = r'(\w+)\s*=\s*([-+]?\d+(?:\.\d+)?|true|false)'
        matches = re.findall(var_pattern, reset_content)
        all_reset_vars.update({var for var, value in matches})

    return all_reset_vars

def extract_header_variables(h_file):
    """从头文件中提取所有类成员变量声明"""
    with open(h_file, 'r', encoding='utf-8') as f:
        content = f.read()
    
    class_pattern = r'class MultiLepPAT\s*:\s*public'
    class_match = re.search(class_pattern, content)
    if not class_match:
        return set()
    
    class_start = class_match.start()
    
    private_start = content.find('private:', class_start)
    if private_start == -1:
        return set()
    
    brace_count = 0
    private_end = private_start
    for i in range(class_start, len(content)):
        if content[i] == '{':
            brace_count += 1
        elif content[i] == '}':
            brace_count -= 1
            if brace_count == 0:
                private_end = i
                break
    
    private_content = content[private_start:private_end]
    
    variable_pattern = r'(?:unsigned\s+int|float|int|bool)\s+(\w+)'
    matches = re.findall(variable_pattern, private_content)
    
    exclude_vars = {
        'isAllVariablesFinite',
    }
    
    return {var for var in matches if var not in exclude_vars}

def main():
    base_dir = Path('/afs/cern.ch/user/q/qinju/cernbox/Onia2MuMu')
    cc_file = base_dir / 'src' / 'MultiLepPAT.cc'
    h_file = base_dir / 'interface' / 'MultiLepPAT.h'
    
    branches = extract_branch_variables(cc_file)
    resets = extract_reset_variables(cc_file)
    headers = extract_header_variables(h_file)
    
    print("=" * 80)
    print("TTree Branch vs resetEventVariables/resetCandidateVariables vs Header 一致性检查")
    print("=" * 80)

    print(f"\n【统计】")
    print(f"  Header 中声明的变量数量: {len(headers)}")
    print(f"  Branch 中定义的变量数量: {len(branches)}")
    print(f"  reset (Event+Candidate) 中重置的变量数量: {len(resets)}")

    print("\n" + "=" * 80)
    print("【问题1】Branch 已定义但 reset 中未重置的变量")
    print("=" * 80)
    branch_not_reset = branches - resets
    if branch_not_reset:
        for var in sorted(branch_not_reset):
            print(f"  ⚠️  {var}")
        print(f"  \n  警告: 以上 {len(branch_not_reset)} 个变量在 Branch 中定义但未在 reset 中重置")
        print(f"        可能导致事件间变量值残留！")
    else:
        print("  ✅ 所有 Branch 变量都已在 reset 中重置")

    print("\n" + "=" * 80)
    print("【问题2】reset 中已重置但未创建 Branch 的变量")
    print("=" * 80)
    reset_not_branch = resets - branches
    if reset_not_branch:
        for var in sorted(reset_not_branch):
            print(f"  ⚠️  {var}")
        print(f"  \n  说明: 以上 {len(reset_not_branch)} 个变量被重置但未写入 TTree")
        print(f"        这些可能是内部计算变量(如 px, py)，不需要保存")
    else:
        print("  ✅ 所有 reset 变量都有对应的 Branch")

    print("\n" + "=" * 80)
    print("【问题3】Header 中声明但既无 Branch 也未 reset 的变量")
    print("=" * 80)
    header_not_used = headers - branches - resets
    if header_not_used:
        for var in sorted(header_not_used):
            print(f"  ⚠️  {var}")
        print(f"  \n  警告: 以上 {len(header_not_used)} 个变量既没有 Branch 也没有 reset")
        print(f"        可能是被遗忘的变量或需要特殊处理的标志位")
    else:
        print("  ✅ 所有 Header 变量都有 Branch 或 reset")

    print("\n" + "=" * 80)
    print("【问题4】Branch 和 reset 中存在但 Header 中未声明的变量")
    print("=" * 80)
    unknown_vars = (branches | resets) - headers
    if unknown_vars:
        for var in sorted(unknown_vars):
            in_branch = "✓ Branch" if var in branches else "✗ Branch"
            in_reset = "✓ reset" if var in resets else "✗ reset"
            print(f"  ⚠️  {var} ({in_branch}, {in_reset})")
        print(f"  \n  错误: 以上 {len(unknown_vars)} 个变量未在 Header 中声明")
        print(f"        可能导致编译错误或运行时问题！")
    else:
        print("  ✅ 所有 Branch 和 reset 变量都在 Header 中声明")
    
    print("\n" + "=" * 80)
    print("【总结】")
    print("=" * 80)
    total_issues = len(branch_not_reset) + len(reset_not_branch) + len(header_not_used) + len(unknown_vars)
    if total_issues == 0:
        print("  🎉 代码一致性检查通过，没有发现问题！")
    else:
        print(f"  发现 {total_issues} 个潜在问题，建议逐一检查")
        if len(branch_not_reset) > 0:
            print(f"    - {len(branch_not_reset)} 个 Branch 变量未在 resetVariables 中重置 [严重]")
        if len(unknown_vars) > 0:
            print(f"    - {len(unknown_vars)} 个变量未在 Header 中声明 [严重]")
        if len(reset_not_branch) > 0:
            print(f"    - {len(reset_not_branch)} 个变量被重置但无 Branch [一般]")
        if len(header_not_used) > 0:
            print(f"    - {len(header_not_used)} 个 Header 变量未使用 [一般]")

if __name__ == '__main__':
    main()