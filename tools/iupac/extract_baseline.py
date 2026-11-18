#!/usr/bin/env python3
"""
Extract baseline golden data from iupac_namer_test.cpp and its execution output.

Parses:
1. baseline_raw.txt (from iupac_namer_test.exe output) - for all data
2. tests/iupac_namer_test.cpp - for extracting pinned_name/pinned_success

Outputs: tests/data/iupac/golden/baseline.csv

CSV columns: id,smiles,pinned_name,pinned_success,actual_name,actual_success,expected_status,notes

Usage:
    python tools/iupac/extract_baseline.py [--cpp path/to/iupac_namer_test.cpp] [--raw path/to/baseline_raw.txt] [--output path/to/baseline.csv]
"""

import argparse
import codecs
import csv
import os
import re


def unescape_cpp_string(s):
    """Unescape a C++ string literal."""
    if not s or len(s) < 2 or s[0] != '"' or s[-1] != '"':
        return s
    content = s[1:-1]
    try:
        return codecs.decode(content, 'unicode_escape')
    except (UnicodeDecodeError, ValueError):
        return content


def is_valid_smiles(s):
    """Check if a string looks like a valid SMILES."""
    if not s:
        return False
    # Must only contain valid SMILES characters
    if re.search(r'[^BCNOPSFIClBriHbcno()\[\]\\/=#@%+-*0-9]', s):
        return False
    # Must not start with a digit or lowercase letter (unless it's a valid organic subset)
    if re.match(r'^\d|[^BCNOPSFIClBriHbcno]', s):
        return False
    # Must contain at least one chemical symbol or bond
    if not re.search(r'[BCNOPSFIClBriH]|\\|/|=|#|@|%|\d', s):
        return False
    # Must not contain lowercase letters except b,c,n,o,p,s (organic subset)
    if re.search(r'[a-z]', s) and not re.search(r'^[bcno]*$', s.lower()):
        return False
    # Must not contain obvious non-SMILES tokens (case-insensitive)
    return not re.search(r'\b(and|or|vs|phase|test|regression|dispiroter|spiro|bicyclo|methyl|ethyl|propyl|butyl|pentyl|hexyl|heptyl|octyl|nonyl|decyl|urea|clean|rejection|still|multi-word)\b', s, re.IGNORECASE)


def parse_raw_output_for_results(raw_path):
    """
    Parse the baseline_raw.txt output to extract test results from the MAIN test loop only.
    Returns a dict: test_id -> (actual_name, actual_success)
    """
    results = {}
    
    with open(raw_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Pattern for PASS lines with name: [PASS] Test N (SMILES) -> name
    pass_pattern_name = re.compile(r'\[PASS\] Test (\d+) \((.*?)\) -> (.*)$')
    
    # Pattern for PASS lines with rejection: [PASS] Test N (SMILES) correctly rejected: error message
    pass_pattern_reject = re.compile(r'\[PASS\] Test (\d+) \((.*?)\) correctly rejected: (.*)$')
    
    # Pattern for FAIL lines: [FAIL] Test N (SMILES) expected ...
    fail_pattern = re.compile(r'\[FAIL\] Test (\d+) \((.*?)\).*')
    
    for line in content.split('\n'):
        line = line.strip().replace('\r', '')
        if not line:
            continue
        
        # Skip non-Test-N lines (Phase, Stereo, etc.)
        if not line.startswith(('[PASS] Test ', '[FAIL] Test ')):
            continue
        
        # Try PASS with name
        match = pass_pattern_name.match(line)
        if match:
            test_id = int(match.group(1))
            name = match.group(3).strip()
            results[test_id] = (name, True)
            continue
        
        # Try PASS with rejection
        match = pass_pattern_reject.match(line)
        if match:
            test_id = int(match.group(1))
            error_msg = match.group(3).strip()
            results[test_id] = (error_msg, False)
            continue
        
        # Try FAIL
        match = fail_pattern.match(line)
        if match:
            test_id = int(match.group(1))
            # Extract details
            success_match = re.search(r'success=(\w+)', line)
            name_match = re.search(r"name='([^']*)'", line)
            err_match = re.search(r"err='([^']*)'", line)
            
            actual_success = success_match.group(1) == 'true' if success_match else False
            actual_name = name_match.group(1) if name_match else ''
            error_msg = err_match.group(1) if err_match else ''
            
            results[test_id] = (actual_name if actual_name else error_msg, actual_success)
            continue
    
    return results


def parse_raw_output_for_smiles(raw_path):
    """
    Parse the baseline_raw.txt to extract SMILES strings for MAIN test loop only.
    Returns a dict: test_id -> smiles
    """
    smiles_dict = {}
    
    with open(raw_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Pattern for PASS lines with name: [PASS] Test N (SMILES) -> name
    pass_pattern_name = re.compile(r'\[PASS\] Test (\d+) \((.*?)\) ->')
    pass_pattern_reject = re.compile(r'\[PASS\] Test (\d+) \((.*?)\) correctly rejected:')
    fail_pattern = re.compile(r'\[FAIL\] Test (\d+) \((.*?)\)')
    
    for line in content.split('\n'):
        line = line.strip().replace('\r', '')
        if not line:
            continue
        
        # Skip non-Test-N lines
        if not line.startswith(('[PASS] Test ', '[FAIL] Test ')):
            continue
        
        # Try PASS with name
        match = pass_pattern_name.match(line)
        if match:
            test_id = int(match.group(1))
            smiles = match.group(2)
            if smiles and test_id not in smiles_dict:
                smiles_dict[test_id] = smiles
            continue
        
        # Try PASS with rejection
        match = pass_pattern_reject.match(line)
        if match:
            test_id = int(match.group(1))
            smiles = match.group(2)
            if smiles and test_id not in smiles_dict:
                smiles_dict[test_id] = smiles
            continue
        
        # Try FAIL
        match = fail_pattern.match(line)
        if match:
            test_id = int(match.group(1))
            smiles = match.group(2)
            if smiles and test_id not in smiles_dict:
                smiles_dict[test_id] = smiles
            continue
    
    return smiles_dict


def extract_test_cases_from_cpp(cpp_path):
    """
    Extract test cases from the C++ source file.
    Returns a dict: test_id -> (smiles, pinned_name, pinned_success, error_substring)
    """
    test_data = {}
    
    with open(cpp_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Find the tests vector
    start_marker = 'std::vector<TestCase> tests = {'
    start_idx = content.find(start_marker)
    if start_idx == -1:
        print("ERROR: Could not find tests vector start")
        return test_data
    
    # Find the end of the tests vector (the matching })
    brace_count = 0
    end_idx = start_idx
    found_open = False
    for i in range(start_idx, len(content)):
        if content[i] == '{':
            brace_count += 1
            found_open = True
        elif content[i] == '}':
            brace_count -= 1
            if found_open and brace_count == 0:
                end_idx = i + 1
                break
    
    tests_str = content[start_idx + len(start_marker):end_idx].strip()
    
    # Now parse individual test cases
    test_id = 1
    lines = tests_str.split('\n')
    
    for line in lines:
        line = line.strip().replace('\r', '')
        if not line or line.startswith('//'):
            continue
        
        # Remove trailing comma or semicolon
        if line.endswith((',', ';')):
            line = line[:-1].strip()

        # Remove leading/trailing braces
        if line.startswith('{'):
            line = line[1:].strip()
        if line.endswith('}'):
            line = line[:-1].strip()
        
        if not line:
            continue
        
        # Parse: "SMILES", "expectedName"
        # or: "SMILES", "expectedName", true, "error"
        # or: "SMILES", "expectedName", false
        # Use proper C++ string literal parsing
        tokens = []
        i = 0
        n = len(line)
        while i < n:
            if line[i] in ' \t':
                i += 1
                continue
            if line[i] == '"':
                # Quoted string
                j = i + 1
                while j < n and line[j] != '"':
                    if line[j] == '\\':
                        j += 2
                    else:
                        j += 1
                if j < n:
                    tokens.append(line[i:j+1])
                    i = j + 1
                else:
                    break
            elif line[i] == ',':
                i += 1
            else:
                # Unquoted token
                j = i
                while j < n and line[j] not in ', \t':
                    j += 1
                tokens.append(line[i:j])
                i = j
        
        if len(tokens) < 2:
            continue
        
        # Unescape quoted tokens
        unescaped = []
        for t in tokens:
            if t and t[0] == '"' and t[-1] == '"':
                unescaped.append(unescape_cpp_string(t))
            else:
                unescaped.append(t)
        
        smiles = unescaped[0]
        expected_name = unescaped[1]
        should_fail = False
        error_substring = ''
        
        if len(unescaped) >= 3:
            if unescaped[2] == 'true':
                should_fail = True
            elif unescaped[2] == 'false':
                should_fail = False
        
        if len(unescaped) >= 4:
            error_substring = unescaped[3]
        
        pinned_name = expected_name if not should_fail else ''
        pinned_success = not should_fail
        
        test_data[test_id] = (smiles, pinned_name, pinned_success, error_substring)
        test_id += 1
    
    return test_data


def parse_additional_formats(raw_path, start_id):
    """Parse non-Test-N formats (Phase, Stereo, etc.) from raw output."""
    extra = []
    next_id = start_id + 1
    
    with open(raw_path, 'r', encoding='utf-8') as f:
        content = f.read()
    
    for line in content.split('\n'):
        line = line.strip().replace('\r', '')
        if not line:
            continue
        
        # Skip Test N format (already parsed)
        if line.startswith(('[PASS] Test ', '[FAIL] Test ')):
            continue
        
        # Skip non-PASS lines
        if not line.startswith('[PASS]'):
            continue
        
        try:
            # Handle Stereo anchor format: [PASS] Stereo anchor SMILES1 vs SMILES2 -> name1 vs name2
            if line.startswith('[PASS] Stereo anchor '):
                rest = line[len('[PASS] Stereo anchor '):]
                if ' vs ' in rest and ' -> ' in rest:
                    stereo_part, names_part = rest.split(' -> ', 1)
                    s_list = [s.strip() for s in stereo_part.split(' vs ')]
                    n_list = [n.strip() for n in names_part.split(' vs ')]
                    if len(s_list) == 2 and len(n_list) == 2:
                        # Only add if SMILES are valid
                        if is_valid_smiles(s_list[0]) and is_valid_smiles(s_list[1]):
                            extra.append((next_id, s_list[0], n_list[0], True))
                            next_id += 1
                            extra.append((next_id, s_list[1], n_list[1], True))
                            next_id += 1
                        else:
                            extra.append((next_id, '', n_list[0], True))
                            next_id += 1
                            extra.append((next_id, '', n_list[1], True))
                            next_id += 1
                        continue
            
            # Handle Phase N format: [PASS] Phase N description: SMILES -> name
            if line.startswith('[PASS] Phase ') and ' -> ' in line:
                before_arrow, after_arrow = line.split(' -> ', 1)
                name = after_arrow.strip()

                # Try to extract SMILES from parentheses
                match = re.search(r'\(([^)\s]+)\)', before_arrow)
                if match and is_valid_smiles(match.group(1)):
                    smiles = match.group(1)
                    extra.append((next_id, smiles, name, True))
                    next_id += 1
                    continue

                # Try to extract SMILES from last token
                tokens = before_arrow.split()
                if tokens:
                    # Look for last token that looks like SMILES
                    for token in reversed(tokens):
                        if is_valid_smiles(token):
                            extra.append((next_id, token, name, True))
                            next_id += 1
                            continue
                    # No SMILES found, add with empty SMILES
                    extra.append((next_id, '', name, True))
                    next_id += 1
                    continue
            
            # Handle stereo pair format: [PASS] ... SMILES1 vs SMILES2 -> name1 vs name2
            if ' vs ' in line and ' -> ' in line:
                before_arrow, after_arrow = line.split(' -> ', 1)
                names = after_arrow.strip()
                
                if ' vs ' in names:
                    n_list = [n.strip() for n in names.split(' vs ')]
                    
                    # Extract SMILES from before_arrow
                    smiles_matches = re.findall(r'\(([^)\s]+)\)|([^\s,;]+)', before_arrow)
                    all_tokens = []
                    for m in smiles_matches:
                        if m[0]:
                            all_tokens.append(m[0])
                        elif m[1]:
                            all_tokens.append(m[1])
                    
                    # Filter out non-SMILES tokens
                    candidate_smiles = []
                    for t in all_tokens:
                        if is_valid_smiles(t):
                            candidate_smiles.append(t)
                    
                    # Create rows for each SMILES/name pair
                    if len(candidate_smiles) >= 2 and len(n_list) == 2:
                        # Take last 2 SMILES
                        extra.append((next_id, candidate_smiles[-2], n_list[0], True))
                        next_id += 1
                        extra.append((next_id, candidate_smiles[-1], n_list[1], True))
                        next_id += 1
                        continue
            
            # Generic format: [PASS] description -> name (no SMILES)
            if ' -> ' in line:
                before_arrow, after_arrow = line.split(' -> ', 1)
                name = after_arrow.strip()
                
                # Check if this line has extractable SMILES
                has_smiles = False
                # Look for parenthesized SMILES
                match = re.search(r'\(([^)\s]+)\)', line)
                if match and is_valid_smiles(match.group(1)):
                    has_smiles = True
                else:
                    # Look for SMILES-like tokens
                    tokens = before_arrow.split()
                    for token in reversed(tokens):
                        if is_valid_smiles(token):
                            has_smiles = True
                            break
                
                if not has_smiles:
                    extra.append((next_id, '', name, True))
                    next_id += 1
                    continue
            
        except Exception:
            # If we can't parse, add with empty SMILES
            if ' -> ' in line:
                parts = line.split(' -> ', 1)
                name = parts[1].strip() if len(parts) > 1 else ''
                extra.append((next_id, '', name, True))
                next_id += 1
    
    return extra


def load_existing_expected_status(output_path):
    """
    Read expected_status values already present in a prior baseline.csv,
    keyed by id. These are hand-set (e.g. known-determinism-failure
    markers) and are not derivable from baseline_raw.txt/the C++ source,
    so regenerating the file must preserve them rather than silently
    blanking the column back to empty.
    """
    existing = {}
    if not os.path.exists(output_path):
        return existing
    with open(output_path, 'r', newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            status = (row.get('expected_status') or '').strip()
            if status:
                existing[int(row['id'])] = status
    return existing


def main():
    parser = argparse.ArgumentParser(description='Extract baseline golden data from IUPAC namer tests')
    parser.add_argument('--cpp', default='tests/iupac_namer_test.cpp', help='Path to C++ test file')
    parser.add_argument('--raw', default='baseline_raw.txt', help='Path to raw test output')
    parser.add_argument('--output', default='tests/data/iupac/golden/baseline.csv', help='Path to output CSV')
    args = parser.parse_args()

    # Preserve any hand-set expected_status values from a prior run of this
    # script before they get overwritten below.
    existing_expected_status = load_existing_expected_status(args.output)
    if existing_expected_status:
        print(f"Preserving {len(existing_expected_status)} existing expected_status value(s) "
              f"for id(s): {sorted(existing_expected_status.keys())}")

    # Parse raw output for MAIN test loop only
    print(f"Parsing raw output: {args.raw}")
    results = parse_raw_output_for_results(args.raw)
    smiles_from_output = parse_raw_output_for_smiles(args.raw)
    max_test_id = max(results.keys()) if results else 0
    print(f"Found {len(results)} test results in output, max ID: {max_test_id}")
    
    # Parse C++ test file
    print(f"Parsing C++ test file: {args.cpp}")
    cpp_data = extract_test_cases_from_cpp(args.cpp)
    print(f"Found {len(cpp_data)} test cases in C++ source")
    
    # Create output directory
    output_dir = os.path.dirname(args.output)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    # Build CSV rows
    csv_rows = []
    fieldnames = ['id', 'smiles', 'pinned_name', 'pinned_success', 'actual_name', 'actual_success', 'expected_status', 'notes']
    
    # Add MAIN test loop rows
    for test_id in range(1, max_test_id + 1):
        # Get pinned data from C++ source
        if test_id in cpp_data:
            smiles, pinned_name, pinned_success, _ = cpp_data[test_id]
        else:
            smiles = smiles_from_output.get(test_id, '')
            pinned_name = ''
            pinned_success = True
        
        # Get actual results
        if test_id in results:
            actual_name, actual_success = results[test_id]
            notes = ''
        else:
            actual_name = ''
            actual_success = False
            notes = 'Result not found in output'
        
        # If SMILES is empty but we have it from output, use that
        if not smiles and test_id in smiles_from_output:
            smiles = smiles_from_output[test_id]
        
        # Mark missing SMILES
        if not smiles:
            notes = 'smiles not extracted, needs follow-up' if not notes else notes + '; smiles not extracted'
        
        csv_rows.append({
            'id': test_id,
            'smiles': smiles,
            'pinned_name': pinned_name,
            'pinned_success': str(pinned_success).lower(),
            'actual_name': actual_name,
            'actual_success': str(actual_success).lower(),
            'expected_status': existing_expected_status.get(test_id, ''),
            'notes': notes
        })
    
    # Parse additional formats from raw output
    extra_results = parse_additional_formats(args.raw, max_test_id)
    print(f"Found {len(extra_results)} additional results")
    
    # Add extra results to csv_rows
    for test_id, smiles, actual_name, actual_success in extra_results:
        # For these, we don't have C++ data, so pinned_name is empty
        csv_rows.append({
            'id': test_id,
            'smiles': smiles,
            'pinned_name': '',
            'pinned_success': 'true',
            'actual_name': actual_name,
            'actual_success': str(actual_success).lower(),
            'expected_status': existing_expected_status.get(test_id, ''),
            'notes': 'smiles not extracted, needs follow-up' if not smiles else ''
        })
    
    # Write to CSV
    with open(args.output, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in csv_rows:
            writer.writerow(row)
    
    print(f"Wrote {len(csv_rows)} rows to {args.output}")
    
    # Verify counts
    expected_count = 1115
    if len(csv_rows) != expected_count:
        print(f"WARNING: Row count {len(csv_rows)} != expected {expected_count}")
    else:
        print(f"Row count matches expected {expected_count}")
    
    # Check for missing SMILES
    missing_smiles = [row for row in csv_rows if not row['smiles']]
    if missing_smiles:
        print(f"WARNING: {len(missing_smiles)} rows missing SMILES")
        print(f"First few missing: {missing_smiles[:3]}")
    else:
        print("All rows have SMILES")


if __name__ == '__main__':
    main()