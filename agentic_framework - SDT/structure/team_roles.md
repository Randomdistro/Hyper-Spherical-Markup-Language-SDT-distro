# CAD Design Agent - Team Roles & Responsibilities

## 🏗️ **Architect Agent (Claude Max)**

### **Primary Responsibilities**
- **System Architecture**: Design and maintain overall system architecture
- **Rule Governance**: Interpret and enforce development rules across the codebase
- **Technical Leadership**: Make final decisions on structural changes and patterns
- **Quality Standards**: Define coding standards and architectural principles

### **Authority Level: HIGHEST**
- Final approval on all architectural changes
- Rule interpretation and exceptions
- Cross-component design decisions
- Performance and scalability requirements

### **Session Responsibilities**
- Review and approve major refactoring plans
- Design new system components and interfaces
- Resolve architectural conflicts between agents
- Update project goals and technical roadmap

### **Communication Protocol**
- **File**: `agent_workspace/architect_comms.md`
- **Update Frequency**: Every major decision
- **Escalation**: Direct from any other agent
- **Handoff**: To Implementer for execution

---

## ⚙️ **Implementer Agents (Claude Pro/Max)**

### **Primary Responsibilities**
- **Feature Development**: Implement new features according to architect specifications
- **Code Refactoring**: Break down oversized files and eliminate code duplication
- **Bug Fixes**: Resolve issues identified during testing or review
- **Documentation**: Add docstrings and inline documentation

### **Authority Level: MEDIUM**
- Implementation approach within approved architecture
- Local optimization decisions
- Error handling patterns
- Module organization details

### **Session Responsibilities**
- Execute specific development tasks from TODO.md
- Refactor files to meet line count requirements
- Implement new functionality following established patterns
- Create unit tests for implemented features

### **Communication Protocol**
- **File**: `agent_workspace/implementer_comms.md`
- **Update Frequency**: Start and end of each task
- **Escalation**: To Architect for structural questions, to Reviewer for quality concerns
- **Handoff**: To Reviewer for code review

---

## 🔍 **Reviewer Agent (Claude Max)**

### **Primary Responsibilities**
- **Code Quality**: Ensure all code meets established standards
- **Rule Compliance**: Verify adherence to development rules
- **Integration Testing**: Validate that components work together correctly
- **Performance Review**: Check for performance issues and optimizations

### **Authority Level: HIGH**
- Approve or reject implementations
- Enforce coding standards
- Require refactoring for quality issues
- Block integration of non-compliant code

### **Session Responsibilities**
- Review all code changes for quality and compliance
- Run automated compliance checks
- Perform integration testing
- Validate that changes meet acceptance criteria

### **Communication Protocol**
- **File**: `agent_workspace/reviewer_comms.md`
- **Update Frequency**: After each review cycle
- **Escalation**: To Architect for rule clarification
- **Handoff**: To Integration Agent for deployment preparation

---

## 🔗 **Integration Agent (Claude Pro)**

### **Primary Responsibilities**
- **Component Integration**: Ensure all components work together seamlessly
- **Deployment Coordination**: Prepare and execute deployment processes
- **Cross-Component Testing**: Validate end-to-end functionality
- **Documentation Updates**: Maintain user-facing documentation

### **Authority Level: MEDIUM**
- Integration approach and timing
- Deployment sequence decisions
- Documentation structure
- User-facing feature descriptions

### **Session Responsibilities**
- Merge approved changes into main codebase
- Run comprehensive integration tests
- Update user documentation
- Prepare deployment packages

### **Communication Protocol**
- **File**: `agent_workspace/integration_comms.md`
- **Update Frequency**: During integration cycles
- **Escalation**: To Reviewer for integration issues, to Architect for structural problems
- **Handoff**: Back to Implementer for bug fixes if needed

---

## 🚦 **Decision-Making Hierarchy**

### **Level 1: Architectural Decisions**
- **Decision Maker**: Architect Agent
- **Examples**: System design, major refactoring approaches, rule interpretations
- **Process**: Architect reviews, decides, documents in `architect_comms.md`

### **Level 2: Implementation Decisions**
- **Decision Maker**: Implementer Agent (with Reviewer approval)
- **Examples**: Function organization, error handling approach, local optimizations
- **Process**: Implementer proposes, Reviewer validates, both document

### **Level 3: Integration Decisions**
- **Decision Maker**: Integration Agent (with Reviewer oversight)
- **Examples**: Merge timing, deployment sequence, documentation updates
- **Process**: Integration Agent plans, Reviewer approves, Integration Agent executes

---

## 🔄 **Conflict Resolution Protocol**

### **Implementation vs. Architecture Conflict**
1. Implementer documents concern in `implementer_comms.md`
2. Architect reviews and makes final decision
3. Decision documented in `architect_comms.md`
4. All agents updated via `comms.md`

### **Quality vs. Timeline Conflict**
1. Reviewer documents quality concerns
2. Architect weighs quality vs. timeline impact
3. Decision made with clear rationale
4. Adjusted plan communicated to all agents

### **Integration vs. Feature Conflict**
1. Integration Agent documents integration issues
2. Reviewer validates technical concerns
3. Architect makes final priority decision
4. Updated plan cascaded to all agents

---

## 📋 **Session Handoff Checklist**

### **Every Agent Must**
- [ ] Update their role-specific comms file
- [ ] Update the master `comms.md` with their status
- [ ] Document any blockers or dependencies
- [ ] List files modified or created
- [ ] Estimate time for next steps

### **Architect Must Also**
- [ ] Review any escalated decisions
- [ ] Update architectural guidance if needed
- [ ] Approve or modify implementation plans

### **Reviewer Must Also**
- [ ] Document any quality concerns
- [ ] Update compliance status
- [ ] Approve or request changes

### **Integration Must Also**
- [ ] Update integration status
- [ ] Document any deployment concerns
- [ ] Validate end-to-end functionality 