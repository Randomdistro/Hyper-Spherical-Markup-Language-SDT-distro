/**
 * HSML Multi-Language Parser
 * AST parsing for HSML, CSSS, ShapeScript, and StyleBot languages
 */
import { TokenType } from './hsml-lexer.js';
// AST Node Types
export var ASTNodeType;
(function (ASTNodeType) {
    // Program structure
    ASTNodeType["PROGRAM"] = "PROGRAM";
    ASTNodeType["MODULE"] = "MODULE";
    ASTNodeType["IMPORT"] = "IMPORT";
    ASTNodeType["EXPORT"] = "EXPORT";
    // HSML-specific nodes
    ASTNodeType["HSML_ELEMENT"] = "HSML_ELEMENT";
    ASTNodeType["HSML_ATTRIBUTE"] = "HSML_ATTRIBUTE";
    ASTNodeType["HSML_CONTENT"] = "HSML_CONTENT";
    // CSSS-specific nodes
    ASTNodeType["CSSS_RULE"] = "CSSS_RULE";
    ASTNodeType["CSSS_SELECTOR"] = "CSSS_SELECTOR";
    ASTNodeType["CSSS_DECLARATION"] = "CSSS_DECLARATION";
    ASTNodeType["CSSS_MATERIAL"] = "CSSS_MATERIAL";
    ASTNodeType["CSSS_ANIMATION"] = "CSSS_ANIMATION";
    ASTNodeType["CSSS_KEYFRAME"] = "CSSS_KEYFRAME";
    // ShapeScript-specific nodes
    ASTNodeType["SHAPE_BEHAVIOR"] = "SHAPE_BEHAVIOR";
    ASTNodeType["SHAPE_PHYSICS"] = "SHAPE_PHYSICS";
    ASTNodeType["SHAPE_FORCE"] = "SHAPE_FORCE";
    ASTNodeType["SHAPE_CONSTRAINT"] = "SHAPE_CONSTRAINT";
    ASTNodeType["SHAPE_EVENT"] = "SHAPE_EVENT";
    // StyleBot-specific nodes
    ASTNodeType["STYB_BOT"] = "STYB_BOT";
    ASTNodeType["STYB_AGENT"] = "STYB_AGENT";
    ASTNodeType["STYB_PARALLEL"] = "STYB_PARALLEL";
    ASTNodeType["STYB_RENDER"] = "STYB_RENDER";
    // Common nodes
    ASTNodeType["EXPRESSION"] = "EXPRESSION";
    ASTNodeType["STATEMENT"] = "STATEMENT";
    ASTNodeType["FUNCTION"] = "FUNCTION";
    ASTNodeType["VARIABLE"] = "VARIABLE";
    ASTNodeType["ASSIGNMENT"] = "ASSIGNMENT";
    ASTNodeType["BINARY_OPERATION"] = "BINARY_OPERATION";
    ASTNodeType["UNARY_OPERATION"] = "UNARY_OPERATION";
    ASTNodeType["CALL_EXPRESSION"] = "CALL_EXPRESSION";
    ASTNodeType["MEMBER_EXPRESSION"] = "MEMBER_EXPRESSION";
    ASTNodeType["LITERAL"] = "LITERAL";
    ASTNodeType["IDENTIFIER"] = "IDENTIFIER";
    // Control flow
    ASTNodeType["IF_STATEMENT"] = "IF_STATEMENT";
    ASTNodeType["WHILE_STATEMENT"] = "WHILE_STATEMENT";
    ASTNodeType["FOR_STATEMENT"] = "FOR_STATEMENT";
    ASTNodeType["RETURN_STATEMENT"] = "RETURN_STATEMENT";
    // Spherical coordinate nodes
    ASTNodeType["SPHERICAL_COORDINATE"] = "SPHERICAL_COORDINATE";
    ASTNodeType["SOLID_ANGLE"] = "SOLID_ANGLE";
    ASTNodeType["MATTER_STATE"] = "MATTER_STATE";
})(ASTNodeType || (ASTNodeType = {}));
// Parser class
export class HSMLMultiParser {
    tokens;
    current = 0;
    language;
    constructor(tokens, language) {
        this.tokens = tokens;
        this.language = language;
    }
    // === MAIN PARSING METHODS ===
    parse() {
        const body = [];
        const imports = [];
        const exports = [];
        while (!this.isAtEnd()) {
            try {
                const node = this.parseTopLevel();
                if (node) {
                    if (node.type === ASTNodeType.IMPORT) {
                        imports.push(node);
                    }
                    else if (node.type === ASTNodeType.EXPORT) {
                        exports.push(node);
                    }
                    else {
                        body.push(node);
                    }
                }
            }
            catch (error) {
                this.synchronize();
            }
        }
        return {
            type: ASTNodeType.PROGRAM,
            body,
            imports,
            exports,
            start: 0,
            end: this.tokens[this.tokens.length - 1]?.column || 0,
            language: this.language
        };
    }
    parseTopLevel() {
        const token = this.peek();
        switch (token.type) {
            case TokenType.IDENTIFIER:
                return this.parseDeclaration();
            case TokenType.LEFT_BRACE:
                return this.parseBlock();
            case TokenType.IF:
                return this.parseIfStatement();
            case TokenType.WHILE:
                return this.parseWhileStatement();
            case TokenType.FOR:
                return this.parseForStatement();
            case TokenType.FUNCTION:
                return this.parseFunction();
            case TokenType.RETURN:
                return this.parseReturnStatement();
            default:
                return this.parseExpression();
        }
    }
    // === HSML PARSING ===
    parseHSMLElement() {
        const start = this.current;
        this.advance(); // consume '<'
        const tagName = this.consume(TokenType.IDENTIFIER, "Expected element tag name").value;
        const attributes = [];
        // Parse attributes
        while (!this.check(TokenType.RIGHT_BRACE) && !this.isAtEnd()) {
            if (this.check(TokenType.IDENTIFIER)) {
                attributes.push(this.parseHSMLAttribute());
            }
            else {
                break;
            }
        }
        const selfClosing = this.check(TokenType.SELF_CLOSING);
        if (selfClosing) {
            this.advance();
            this.advance(); // consume '>'
            return {
                type: ASTNodeType.HSML_ELEMENT,
                tagName,
                attributes,
                children: [],
                selfClosing: true,
                start,
                end: this.current,
                language: 'hsml'
            };
        }
        this.advance(); // consume '>'
        const children = [];
        let content;
        // Parse content and children
        while (!this.check(TokenType.CLOSE_TAG) && !this.isAtEnd()) {
            if (this.check(TokenType.LEFT_BRACE)) {
                children.push(this.parseHSMLElement());
            }
            else {
                // Parse text content
                const textToken = this.advance();
                if (textToken.type === TokenType.STRING) {
                    content = textToken.value;
                }
            }
        }
        if (this.check(TokenType.CLOSE_TAG)) {
            this.advance(); // consume '</'
            this.advance(); // consume tag name
            this.advance(); // consume '>'
        }
        return {
            type: ASTNodeType.HSML_ELEMENT,
            tagName,
            attributes,
            children,
            selfClosing: false,
            content,
            start,
            end: this.current,
            language: 'hsml'
        };
    }
    parseHSMLAttribute() {
        const start = this.current;
        const name = this.advance().value;
        this.consume(TokenType.COLON, "Expected ':' after attribute name");
        const value = this.parseExpression();
        return {
            type: ASTNodeType.HSML_ATTRIBUTE,
            name,
            value,
            start,
            end: this.current,
            language: 'hsml'
        };
    }
    // === CSSS PARSING ===
    parseCSSSRule() {
        const start = this.current;
        const selectors = [];
        // Parse selectors
        while (!this.check(TokenType.LEFT_BRACE) && !this.isAtEnd()) {
            selectors.push(this.parseCSSSSelector());
            if (this.check(TokenType.COMMA)) {
                this.advance();
            }
        }
        this.consume(TokenType.LEFT_BRACE, "Expected '{' after selectors");
        const declarations = [];
        while (!this.check(TokenType.RIGHT_BRACE) && !this.isAtEnd()) {
            declarations.push(this.parseCSSSDeclaration());
        }
        this.consume(TokenType.RIGHT_BRACE, "Expected '}' after declarations");
        return {
            type: ASTNodeType.CSSS_RULE,
            selectors,
            declarations,
            start,
            end: this.current,
            language: 'csss'
        };
    }
    parseCSSSSelector() {
        const start = this.current;
        const selector = this.advance().value;
        // Calculate specificity (simplified)
        const specificity = selector.split(/[.#]/).length;
        return {
            type: ASTNodeType.CSSS_SELECTOR,
            selector,
            specificity,
            start,
            end: this.current,
            language: 'csss'
        };
    }
    parseCSSSDeclaration() {
        const start = this.current;
        const property = this.advance().value;
        this.consume(TokenType.COLON, "Expected ':' after property name");
        const value = this.parseExpression();
        let important = false;
        if (this.check(TokenType.IDENTIFIER) && this.peek().value === 'important') {
            this.advance();
            important = true;
        }
        this.consume(TokenType.SEMICOLON, "Expected ';' after declaration");
        return {
            type: ASTNodeType.CSSS_DECLARATION,
            property,
            value,
            important,
            start,
            end: this.current,
            language: 'csss'
        };
    }
    // === ShapeScript PARSING ===
    parseShapeBehavior() {
        const start = this.current;
        this.consume(TokenType.BEHAVIOR, "Expected 'behavior' keyword");
        const name = this.consume(TokenType.IDENTIFIER, "Expected behavior name").value;
        this.consume(TokenType.LEFT_BRACE, "Expected '{' after behavior name");
        const physics = [];
        const events = [];
        while (!this.check(TokenType.RIGHT_BRACE) && !this.isAtEnd()) {
            if (this.check(TokenType.PHYSICS)) {
                physics.push(this.parseShapePhysics());
            }
            else if (this.check(TokenType.EVENT)) {
                events.push(this.parseShapeEvent());
            }
            else {
                this.advance(); // skip unknown tokens
            }
        }
        this.consume(TokenType.RIGHT_BRACE, "Expected '}' after behavior body");
        return {
            type: ASTNodeType.SHAPE_BEHAVIOR,
            name,
            physics,
            events,
            start,
            end: this.current,
            language: 'shape'
        };
    }
    parseShapePhysics() {
        const start = this.current;
        this.consume(TokenType.PHYSICS, "Expected 'physics' keyword");
        this.consume(TokenType.LEFT_BRACE, "Expected '{' after physics");
        const forces = [];
        const constraints = [];
        let matterState;
        while (!this.check(TokenType.RIGHT_BRACE) && !this.isAtEnd()) {
            if (this.check(TokenType.FORCE)) {
                forces.push(this.parseShapeForce());
            }
            else if (this.check(TokenType.CONSTRAINT)) {
                constraints.push(this.parseShapeConstraint());
            }
            else if (this.isMatterState(this.peek().value)) {
                matterState = this.parseMatterState();
            }
            else {
                this.advance();
            }
        }
        this.consume(TokenType.RIGHT_BRACE, "Expected '}' after physics body");
        return {
            type: ASTNodeType.SHAPE_PHYSICS,
            forces,
            constraints,
            matterState: matterState,
            start,
            end: this.current,
            language: 'shape'
        };
    }
    parseShapeForce() {
        const start = this.current;
        this.consume(TokenType.FORCE, "Expected 'force' keyword");
        const forceType = this.consume(TokenType.IDENTIFIER, "Expected force type").value;
        this.consume(TokenType.LEFT_BRACE, "Expected '{' after force type");
        this.consume(TokenType.IDENTIFIER, "Expected 'magnitude'");
        this.consume(TokenType.COLON, "Expected ':' after magnitude");
        const magnitude = this.parseExpression();
        this.consume(TokenType.COMMA, "Expected ',' after magnitude");
        this.consume(TokenType.IDENTIFIER, "Expected 'direction'");
        this.consume(TokenType.COLON, "Expected ':' after direction");
        const direction = this.parseSphericalCoordinate();
        this.consume(TokenType.RIGHT_BRACE, "Expected '}' after force properties");
        return {
            type: ASTNodeType.SHAPE_FORCE,
            forceType,
            magnitude,
            direction,
            start,
            end: this.current,
            language: 'shape'
        };
    }
    parseShapeConstraint() {
        const start = this.current;
        this.consume(TokenType.CONSTRAINT, "Expected 'constraint' keyword");
        const constraintType = this.consume(TokenType.IDENTIFIER, "Expected constraint type").value;
        this.consume(TokenType.LEFT_BRACE, "Expected '{' after constraint type");
        const parameters = {};
        while (!this.check(TokenType.RIGHT_BRACE) && !this.isAtEnd()) {
            const paramName = this.consume(TokenType.IDENTIFIER, "Expected parameter name").value;
            this.consume(TokenType.COLON, "Expected ':' after parameter name");
            const paramValue = this.parseExpression();
            parameters[paramName] = paramValue;
            if (this.check(TokenType.COMMA)) {
                this.advance();
            }
        }
        this.consume(TokenType.RIGHT_BRACE, "Expected '}' after constraint parameters");
        return {
            type: ASTNodeType.SHAPE_CONSTRAINT,
            constraintType,
            parameters,
            start,
            end: this.current,
            language: 'shape'
        };
    }
    parseShapeEvent() {
        const start = this.current;
        this.consume(TokenType.EVENT, "Expected 'event' keyword");
        const trigger = this.consume(TokenType.IDENTIFIER, "Expected event trigger").value;
        this.consume(TokenType.COLON, "Expected ':' after trigger");
        const response = this.parseExpression();
        return {
            type: ASTNodeType.SHAPE_EVENT,
            trigger,
            response,
            start,
            end: this.current,
            language: 'shape'
        };
    }
    // === StyleBot PARSING ===
    parseStyleBot() {
        const start = this.current;
        this.consume(TokenType.BOT, "Expected 'bot' keyword");
        const name = this.consume(TokenType.IDENTIFIER, "Expected bot name").value;
        this.consume(TokenType.LEFT_BRACE, "Expected '{' after bot name");
        const agents = [];
        let parallel = false;
        while (!this.check(TokenType.RIGHT_BRACE) && !this.isAtEnd()) {
            if (this.check(TokenType.PARALLEL)) {
                this.advance();
                parallel = true;
            }
            else if (this.check(TokenType.AGENT)) {
                agents.push(this.parseStyleBotAgent());
            }
            else {
                this.advance();
            }
        }
        this.consume(TokenType.RIGHT_BRACE, "Expected '}' after bot body");
        return {
            type: ASTNodeType.STYB_BOT,
            name,
            agents,
            parallel,
            start,
            end: this.current,
            language: 'styb'
        };
    }
    parseStyleBotAgent() {
        const start = this.current;
        this.consume(TokenType.AGENT, "Expected 'agent' keyword");
        const name = this.consume(TokenType.IDENTIFIER, "Expected agent name").value;
        this.consume(TokenType.LEFT_BRACE, "Expected '{' after agent name");
        const render = [];
        const optimize = [];
        while (!this.check(TokenType.RIGHT_BRACE) && !this.isAtEnd()) {
            if (this.check(TokenType.RENDER)) {
                render.push(this.parseStyleBotRender());
            }
            else if (this.check(TokenType.OPTIMIZE)) {
                optimize.push(this.parseExpression());
            }
            else {
                this.advance();
            }
        }
        this.consume(TokenType.RIGHT_BRACE, "Expected '}' after agent body");
        return {
            type: ASTNodeType.STYB_AGENT,
            name,
            render,
            optimize,
            start,
            end: this.current,
            language: 'styb'
        };
    }
    parseStyleBotRender() {
        const start = this.current;
        this.consume(TokenType.RENDER, "Expected 'render' keyword");
        const target = this.consume(TokenType.IDENTIFIER, "Expected render target").value;
        this.consume(TokenType.LEFT_BRACE, "Expected '{' after render target");
        this.consume(TokenType.IDENTIFIER, "Expected 'quality'");
        this.consume(TokenType.COLON, "Expected ':' after quality");
        const quality = this.parseExpression();
        this.consume(TokenType.COMMA, "Expected ',' after quality");
        this.consume(TokenType.IDENTIFIER, "Expected 'priority'");
        this.consume(TokenType.COLON, "Expected ':' after priority");
        const priority = parseInt(this.consume(TokenType.NUMBER, "Expected priority number").value);
        this.consume(TokenType.RIGHT_BRACE, "Expected '}' after render properties");
        return {
            type: ASTNodeType.STYB_RENDER,
            target,
            quality,
            priority,
            start,
            end: this.current,
            language: 'styb'
        };
    }
    // === EXPRESSION PARSING ===
    parseExpression() {
        return {
            type: ASTNodeType.EXPRESSION,
            expression: this.parseEquality(),
            start: this.current,
            end: this.current,
            language: this.language
        };
    }
    parseEquality() {
        let expr = this.parseComparison();
        while (this.match(TokenType.EQUALS, TokenType.NOT_EQUALS)) {
            const operator = this.previous().value;
            const right = this.parseComparison();
            expr = {
                type: ASTNodeType.BINARY_OPERATION,
                operator,
                left: { type: ASTNodeType.EXPRESSION, expression: expr, start: expr.start, end: expr.end, language: this.language },
                right: { type: ASTNodeType.EXPRESSION, expression: right, start: right.start, end: right.end, language: this.language },
                start: expr.start,
                end: right.end,
                language: this.language
            };
        }
        return expr;
    }
    parseComparison() {
        let expr = this.parseTerm();
        while (this.match(TokenType.LESS_THAN, TokenType.GREATER_THAN, TokenType.LESS_EQUAL, TokenType.GREATER_EQUAL)) {
            const operator = this.previous().value;
            const right = this.parseTerm();
            expr = {
                type: ASTNodeType.BINARY_OPERATION,
                operator,
                left: { type: ASTNodeType.EXPRESSION, expression: expr, start: expr.start, end: expr.end, language: this.language },
                right: { type: ASTNodeType.EXPRESSION, expression: right, start: right.start, end: right.end, language: this.language },
                start: expr.start,
                end: right.end,
                language: this.language
            };
        }
        return expr;
    }
    parseTerm() {
        let expr = this.parseFactor();
        while (this.match(TokenType.PLUS, TokenType.MINUS)) {
            const operator = this.previous().value;
            const right = this.parseFactor();
            expr = {
                type: ASTNodeType.BINARY_OPERATION,
                operator,
                left: { type: ASTNodeType.EXPRESSION, expression: expr, start: expr.start, end: expr.end, language: this.language },
                right: { type: ASTNodeType.EXPRESSION, expression: right, start: right.start, end: right.end, language: this.language },
                start: expr.start,
                end: right.end,
                language: this.language
            };
        }
        return expr;
    }
    parseFactor() {
        let expr = this.parseUnary();
        while (this.match(TokenType.MULTIPLY, TokenType.DIVIDE)) {
            const operator = this.previous().value;
            const right = this.parseUnary();
            expr = {
                type: ASTNodeType.BINARY_OPERATION,
                operator,
                left: { type: ASTNodeType.EXPRESSION, expression: expr, start: expr.start, end: expr.end, language: this.language },
                right: { type: ASTNodeType.EXPRESSION, expression: right, start: right.start, end: right.end, language: this.language },
                start: expr.start,
                end: right.end,
                language: this.language
            };
        }
        return expr;
    }
    parseUnary() {
        if (this.match(TokenType.MINUS)) {
            const operator = this.previous().value;
            const right = this.parseUnary();
            return {
                type: ASTNodeType.UNARY_OPERATION,
                operator,
                operand: { type: ASTNodeType.EXPRESSION, expression: right, start: right.start, end: right.end, language: this.language },
                start: this.current,
                end: right.end,
                language: this.language
            };
        }
        return this.parseCall();
    }
    parseCall() {
        let expr = this.parsePrimary();
        while (true) {
            if (this.match(TokenType.LEFT_PAREN)) {
                expr = this.finishCall(expr);
            }
            else if (this.match(TokenType.DOT)) {
                const name = this.consume(TokenType.IDENTIFIER, "Expected property name after '.'").value;
                expr = {
                    type: ASTNodeType.MEMBER_EXPRESSION,
                    object: { type: ASTNodeType.EXPRESSION, expression: expr, start: expr.start, end: expr.end, language: this.language },
                    property: name,
                    start: expr.start,
                    end: this.current,
                    language: this.language
                };
            }
            else {
                break;
            }
        }
        return expr;
    }
    finishCall(callee) {
        const arguments_ = [];
        if (!this.check(TokenType.RIGHT_PAREN)) {
            do {
                arguments_.push(this.parseExpression());
            } while (this.match(TokenType.COMMA));
        }
        this.consume(TokenType.RIGHT_PAREN, "Expected ')' after arguments");
        return {
            type: ASTNodeType.CALL_EXPRESSION,
            callee: { type: ASTNodeType.EXPRESSION, expression: callee, start: callee.start, end: callee.end, language: this.language },
            arguments: arguments_,
            start: callee.start,
            end: this.current,
            language: this.language
        };
    }
    parsePrimary() {
        if (this.match(TokenType.NUMBER, TokenType.STRING, TokenType.BOOLEAN)) {
            const token = this.previous();
            return {
                type: ASTNodeType.LITERAL,
                value: token.value,
                raw: token.value,
                start: token.column,
                end: token.column + token.value.length,
                language: this.language
            };
        }
        if (this.match(TokenType.IDENTIFIER)) {
            const token = this.previous();
            return {
                type: ASTNodeType.IDENTIFIER,
                name: token.value,
                start: token.column,
                end: token.column + token.value.length,
                language: this.language
            };
        }
        throw new Error(`Unexpected token: ${this.peek().value}`);
    }
    parseSphericalPrimary() {
        if (this.match(TokenType.SPHERICAL_COORD)) {
            return this.parseSphericalCoordinate();
        }
        if (this.match(TokenType.SOLID_ANGLE)) {
            return this.parseSolidAngle();
        }
        if (this.isMatterState(this.peek().value)) {
            return this.parseMatterState();
        }
        throw new Error(`Unexpected token: ${this.peek().value}`);
    }
    // === SPHERICAL COORDINATE PARSING ===
    parseSphericalCoordinate() {
        const start = this.current;
        this.consume(TokenType.LEFT_PAREN, "Expected '(' for spherical coordinate");
        this.consume(TokenType.IDENTIFIER, "Expected 'r'");
        this.consume(TokenType.COLON, "Expected ':' after 'r'");
        const r = this.parseExpression();
        this.consume(TokenType.COMMA, "Expected ',' after r value");
        this.consume(TokenType.IDENTIFIER, "Expected 'θ' or 'theta'");
        this.consume(TokenType.COLON, "Expected ':' after theta");
        const theta = this.parseExpression();
        this.consume(TokenType.COMMA, "Expected ',' after theta value");
        this.consume(TokenType.IDENTIFIER, "Expected 'φ' or 'phi'");
        this.consume(TokenType.COLON, "Expected ':' after phi");
        const phi = this.parseExpression();
        this.consume(TokenType.RIGHT_PAREN, "Expected ')' for spherical coordinate");
        return {
            type: ASTNodeType.SPHERICAL_COORDINATE,
            r,
            theta,
            phi,
            start,
            end: this.current,
            language: this.language
        };
    }
    parseSolidAngle() {
        const start = this.current;
        this.consume(TokenType.SOLID_ANGLE, "Expected 'Ω' for solid angle");
        this.consume(TokenType.LEFT_PAREN, "Expected '(' for solid angle");
        this.consume(TokenType.IDENTIFIER, "Expected 'ω' or 'omega'");
        this.consume(TokenType.COLON, "Expected ':' after omega");
        const omega = this.parseExpression();
        this.consume(TokenType.COMMA, "Expected ',' after omega value");
        this.consume(TokenType.IDENTIFIER, "Expected 'θ_min'");
        this.consume(TokenType.COLON, "Expected ':' after theta_min");
        const theta_min = this.parseExpression();
        this.consume(TokenType.COMMA, "Expected ',' after theta_min value");
        this.consume(TokenType.IDENTIFIER, "Expected 'θ_max'");
        this.consume(TokenType.COLON, "Expected ':' after theta_max");
        const theta_max = this.parseExpression();
        this.consume(TokenType.COMMA, "Expected ',' after theta_max value");
        this.consume(TokenType.IDENTIFIER, "Expected 'φ_min'");
        this.consume(TokenType.COLON, "Expected ':' after phi_min");
        const phi_min = this.parseExpression();
        this.consume(TokenType.COMMA, "Expected ',' after phi_min value");
        this.consume(TokenType.IDENTIFIER, "Expected 'φ_max'");
        this.consume(TokenType.COLON, "Expected ':' after phi_max");
        const phi_max = this.parseExpression();
        this.consume(TokenType.RIGHT_PAREN, "Expected ')' for solid angle");
        return {
            type: ASTNodeType.SOLID_ANGLE,
            omega,
            theta_min,
            theta_max,
            phi_min,
            phi_max,
            start,
            end: this.current,
            language: this.language
        };
    }
    parseMatterState() {
        const start = this.current;
        const stateToken = this.advance();
        const state = stateToken.value;
        const properties = {};
        if (this.match(TokenType.LEFT_BRACE)) {
            while (!this.check(TokenType.RIGHT_BRACE) && !this.isAtEnd()) {
                const propertyName = this.consume(TokenType.IDENTIFIER, "Expected property name").value;
                this.consume(TokenType.COLON, "Expected ':' after property name");
                const propertyValue = this.parseExpression();
                properties[propertyName] = propertyValue;
                if (this.check(TokenType.COMMA)) {
                    this.advance();
                }
            }
            this.consume(TokenType.RIGHT_BRACE, "Expected '}' after matter state properties");
        }
        return {
            type: ASTNodeType.MATTER_STATE,
            state,
            properties,
            start,
            end: this.current,
            language: this.language
        };
    }
    // === UTILITY METHODS ===
    match(...types) {
        for (const type of types) {
            if (this.check(type)) {
                this.advance();
                return true;
            }
        }
        return false;
    }
    check(type) {
        if (this.isAtEnd())
            return false;
        return this.peek().type === type;
    }
    advance() {
        if (!this.isAtEnd())
            this.current++;
        return this.previous();
    }
    isAtEnd() {
        return this.peek().type === TokenType.EOF;
    }
    peek() {
        return this.tokens[this.current];
    }
    previous() {
        return this.tokens[this.current - 1];
    }
    consume(type, message) {
        if (this.check(type))
            return this.advance();
        throw new Error(message);
    }
    isMatterState(value) {
        return ['solid', 'liquid', 'gas', 'plasma'].includes(value);
    }
    synchronize() {
        this.advance();
        while (!this.isAtEnd()) {
            if (this.previous().type === TokenType.SEMICOLON)
                return;
            switch (this.peek().type) {
                case TokenType.FUNCTION:
                case TokenType.VAR:
                case TokenType.LET:
                case TokenType.CONST:
                case TokenType.IF:
                case TokenType.WHILE:
                case TokenType.FOR:
                case TokenType.RETURN:
                    return;
            }
            this.advance();
        }
    }
    // Additional parsing methods for control flow and declarations
    parseDeclaration() {
        // Implementation for variable and function declarations
        return this.parseExpression();
    }
    parseBlock() {
        // Implementation for block statements
        return this.parseExpression();
    }
    parseIfStatement() {
        // Implementation for if statements
        return {
            type: ASTNodeType.IF_STATEMENT,
            condition: this.parseExpression(),
            consequent: this.parseExpression(),
            start: this.current,
            end: this.current,
            language: this.language
        };
    }
    parseWhileStatement() {
        // Implementation for while statements
        return {
            type: ASTNodeType.WHILE_STATEMENT,
            condition: this.parseExpression(),
            body: this.parseExpression(),
            start: this.current,
            end: this.current,
            language: this.language
        };
    }
    parseForStatement() {
        // Implementation for for statements
        return {
            type: ASTNodeType.FOR_STATEMENT,
            init: undefined,
            condition: undefined,
            update: undefined,
            body: this.parseExpression(),
            start: this.current,
            end: this.current,
            language: this.language
        };
    }
    parseFunction() {
        // Implementation for function declarations
        return this.parseExpression();
    }
    parseReturnStatement() {
        // Implementation for return statements
        return {
            type: ASTNodeType.RETURN_STATEMENT,
            argument: undefined,
            start: this.current,
            end: this.current,
            language: this.language
        };
    }
}
//# sourceMappingURL=hsml-parser.js.map