# HSML-SDT Validation Framework - Quick Start Guide

**Get started with validation testing in 5 minutes!**

---

## What is the Validation Framework?

The HSML-SDT Validation Framework proves the correctness of the Spatial Displacement Theory implementation through comprehensive, empirical testing. It validates:

- **Coordinate System**: Zero-division safety, accuracy, edge cases
- **Physics**: Stability, energy conservation, analytical solutions
- **Performance**: Frame times, scalability, memory usage

---

## Prerequisites

```bash
# Node.js 18+ required
node --version  # Should be >= 18.0.0

# Install dependencies
npm install
```

---

## Running Validation Tests

### Quick Validation (All Tests)

```bash
npm run validate:all
```

This runs all three validation suites:
1. Coordinate system validation
2. Physics accuracy validation
3. Performance benchmarks

**Expected runtime**: ~2-5 minutes

### Individual Validation Suites

**Coordinate System** (~30 seconds):
```bash
npm run validate:coordinates
```

Tests:
- 1,000,000+ coordinate operations for zero-division safety
- Distance calculation accuracy
- Edge case handling (poles, boundaries, NaN/Infinity)
- Trigonometric function validity

**Physics Accuracy** (~2-3 minutes):
```bash
npm run validate:physics
```

Tests:
- 100,000+ frame stability simulations
- Energy conservation (<0.1% drift)
- Analytical solution matching
- Multi-body interactions
- Matter state transitions

**Performance Benchmarks** (~1-2 minutes):
```bash
npm run validate:performance
```

Tests:
- Coordinate operation benchmarks
- 21D state evolution performance
- Multi-body system scaling (10, 100, 1000 entities)
- Memory leak detection

### Watch Mode (Development)

```bash
npm run validate:watch
```

Automatically reruns tests when files change. Great for development!

---

## Understanding Test Results

### Success Output

```
=== COORDINATE SYSTEM VALIDATION REPORT ===
Status: TESTS DEFINED
Zero-division tests: ✅ Implemented
Accuracy tests: ✅ Implemented
Edge case tests: ✅ Implemented
Performance tests: ✅ Implemented
==========================================

 PASS  validation/coordinate-system-validation.test.ts
  ✓ should never produce zero angles (1234ms)
  ✓ should handle extreme angle values (45ms)
  ✓ should prevent division by zero (23ms)
  ...

Tests: 52 passed, 52 total
Time: 12.345s
```

### Failure Output

If tests fail, you'll see detailed information:

```
 FAIL  validation/physics-accuracy-validation.test.ts
  ✕ should conserve energy within 0.1% (2345ms)

  Expected: 0.0009
  Received: 0.0015

  Energy drift 0.15% exceeds tolerance 0.1%

  at Object.<anonymous> (validation/physics-accuracy-validation.test.ts:123:7)
```

---

## Performance Benchmarks

Performance tests output detailed statistics:

```
=== PERFORMANCE BENCHMARK RESULTS ===

Safe Coordinate Conversion:
  Iterations: 100000
  Avg: 0.0004ms
  Min: 0.0002ms
  Max: 0.0012ms
  P95: 0.0006ms
  P99: 0.0008ms
  Ops/sec: 2500000

100-Entity Frame:
  Iterations: 100
  Avg: 4.23ms
  Min: 3.89ms
  Max: 5.67ms
  P95: 4.91ms
  P99: 5.34ms
  Ops/sec: 236
```

**Targets**:
- Frame time: <7ms @ 144Hz ✅
- 21D evolution: <50μs ✅
- Coordinate ops: <1μs ✅

---

## Continuous Integration

Validation tests run automatically on:

- **Every push** to main/develop
- **Every pull request**
- **Nightly** (full validation suite)

### CI/CD Workflow

```yaml
# .github/workflows/ci.yml

jobs:
  validation-tests:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - run: npm ci
      - run: npm run validate:all
```

PRs must pass all validation tests before merging.

---

## Line Limit Checking

Check that all files comply with the 1350-line limit:

```bash
npm run check:lines
```

Output:
```
============================================================
HSML-SDT LINE LIMIT CHECK
============================================================

Line Limit: 1350 lines per file
Files Scanned: 47
Compliant: 47
Violations: 0

────────────────────────────────────────────────────────────
✅ SUCCESS: All files comply with the line limit!
────────────────────────────────────────────────────────────
```

---

## Generating Validation Reports

Create a JSON report for analysis:

```bash
npm run validate:report
```

Outputs to `validation-report.json`:

```json
{
  "testResults": [...],
  "numPassedTests": 152,
  "numFailedTests": 0,
  "success": true,
  "startTime": 1697234567890,
  "endTime": 1697234612345
}
```

---

## Common Issues

### Issue: Tests Timeout

**Symptom**: Tests take too long and timeout

**Solution**: Increase Jest timeout:

```javascript
// In test file
test('long running test', () => {
  // ...
}, 60000); // 60 second timeout
```

### Issue: Performance Tests Fail on Slow Machine

**Symptom**: Performance benchmarks don't meet targets

**Solution**: Performance targets are calibrated for modern hardware. On slower machines, adjust targets in validation tests or skip performance validation:

```bash
npm run validate:coordinates && npm run validate:physics
```

### Issue: Memory Errors

**Symptom**: "JavaScript heap out of memory"

**Solution**: Increase Node memory:

```bash
NODE_OPTIONS=--max-old-space-size=4096 npm run validate:all
```

---

## Advanced Usage

### Custom Test Patterns

Run specific tests:

```bash
# Run only zero-division tests
npm test -- --testNamePattern="zero.division"

# Run tests in a specific file
npm test validation/coordinate-system-validation.test.ts
```

### Debugging Tests

```bash
# Run with verbose output
npm test -- --verbose

# Run in watch mode with coverage
npm test -- --watch --coverage
```

### Benchmarking Specific Operations

```bash
# Run only performance benchmarks
npm run validate:performance

# Compare performance over time
npm run validate:performance > benchmarks-$(date +%Y%m%d).txt
```

---

## Integration with Development Workflow

### Pre-commit Hook (Recommended)

Add to `.git/hooks/pre-commit`:

```bash
#!/bin/sh
npm run check:lines
npm run typecheck
```

### Pre-push Hook

Add to `.git/hooks/pre-push`:

```bash
#!/bin/sh
npm run validate:all
```

---

## Next Steps

1. **Run your first validation**: `npm run validate:all`
2. **Review the results**: Check for any failures
3. **Set up CI/CD**: Push to GitHub to trigger automated testing
4. **Enable watch mode**: `npm run validate:watch` while developing

---

## Help & Support

**Documentation**:
- Full validation strategy: `validation/README.md`
- Test utilities: `validation/utils/validation-helpers.ts`
- Project goals: `PROJECT_GOALS_HSML_SDT.md`

**Framework Rules**:
- Line limits: `agentic_framework - SDT/rules/2_line_limit.md`
- Testing requirements: `agentic_framework - SDT/rules/5_error_limiting.md`

**Questions?**
- Check existing tests for examples
- Review validation helper utilities
- See PROJECT_GOALS_HSML_SDT.md for acceptance criteria

---

**Validation Framework Status**: ✅ READY FOR USE
**Last Updated**: 2025-10-14
**Version**: 1.0.0
