/**
 * HSML Multi-Language Parser
 * AST parsing for HSML, CSSS, ShapeScript, and StyleBot languages
 */
import { Token } from './hsml-lexer.js';
export declare enum ASTNodeType {
    PROGRAM = "PROGRAM",
    MODULE = "MODULE",
    IMPORT = "IMPORT",
    EXPORT = "EXPORT",
    HSML_ELEMENT = "HSML_ELEMENT",
    HSML_ATTRIBUTE = "HSML_ATTRIBUTE",
    HSML_CONTENT = "HSML_CONTENT",
    CSSS_RULE = "CSSS_RULE",
    CSSS_SELECTOR = "CSSS_SELECTOR",
    CSSS_DECLARATION = "CSSS_DECLARATION",
    CSSS_MATERIAL = "CSSS_MATERIAL",
    CSSS_ANIMATION = "CSSS_ANIMATION",
    CSSS_KEYFRAME = "CSSS_KEYFRAME",
    SHAPE_BEHAVIOR = "SHAPE_BEHAVIOR",
    SHAPE_PHYSICS = "SHAPE_PHYSICS",
    SHAPE_FORCE = "SHAPE_FORCE",
    SHAPE_CONSTRAINT = "SHAPE_CONSTRAINT",
    SHAPE_EVENT = "SHAPE_EVENT",
    STYB_BOT = "STYB_BOT",
    STYB_AGENT = "STYB_AGENT",
    STYB_PARALLEL = "STYB_PARALLEL",
    STYB_RENDER = "STYB_RENDER",
    EXPRESSION = "EXPRESSION",
    STATEMENT = "STATEMENT",
    FUNCTION = "FUNCTION",
    VARIABLE = "VARIABLE",
    ASSIGNMENT = "ASSIGNMENT",
    BINARY_OPERATION = "BINARY_OPERATION",
    UNARY_OPERATION = "UNARY_OPERATION",
    CALL_EXPRESSION = "CALL_EXPRESSION",
    MEMBER_EXPRESSION = "MEMBER_EXPRESSION",
    LITERAL = "LITERAL",
    IDENTIFIER = "IDENTIFIER",
    IF_STATEMENT = "IF_STATEMENT",
    WHILE_STATEMENT = "WHILE_STATEMENT",
    FOR_STATEMENT = "FOR_STATEMENT",
    RETURN_STATEMENT = "RETURN_STATEMENT",
    SPHERICAL_COORDINATE = "SPHERICAL_COORDINATE",
    SOLID_ANGLE = "SOLID_ANGLE",
    MATTER_STATE = "MATTER_STATE"
}
export interface ASTNode {
    type: ASTNodeType;
    start: number;
    end: number;
    language: 'hsml' | 'csss' | 'shape' | 'styb';
}
export interface ProgramNode extends ASTNode {
    type: ASTNodeType.PROGRAM;
    body: ASTNode[];
    imports: ImportNode[];
    exports: ExportNode[];
}
export interface ImportNode extends ASTNode {
    type: ASTNodeType.IMPORT;
    module: string;
    imports: string[];
    alias?: string;
}
export interface ExportNode extends ASTNode {
    type: ASTNodeType.EXPORT;
    declaration: ASTNode;
}
export interface HSMLElementNode extends ASTNode {
    type: ASTNodeType.HSML_ELEMENT;
    tagName: string;
    attributes: HSMLAttributeNode[];
    children: HSMLElementNode[];
    selfClosing: boolean;
    content?: string;
}
export interface HSMLAttributeNode extends ASTNode {
    type: ASTNodeType.HSML_ATTRIBUTE;
    name: string;
    value: ExpressionNode;
}
export interface CSSSRuleNode extends ASTNode {
    type: ASTNodeType.CSSS_RULE;
    selectors: CSSSSelectorNode[];
    declarations: CSSSDeclarationNode[];
}
export interface CSSSSelectorNode extends ASTNode {
    type: ASTNodeType.CSSS_SELECTOR;
    selector: string;
    specificity: number;
}
export interface CSSSDeclarationNode extends ASTNode {
    type: ASTNodeType.CSSS_DECLARATION;
    property: string;
    value: ExpressionNode;
    important: boolean;
}
export interface CSSSMaterialNode extends ASTNode {
    type: ASTNodeType.CSSS_MATERIAL;
    name: string;
    properties: CSSSDeclarationNode[];
    matterState: MatterStateNode;
}
export interface CSSSAnimationNode extends ASTNode {
    type: ASTNodeType.CSSS_ANIMATION;
    name: string;
    keyframes: CSSSKeyframeNode[];
    duration: ExpressionNode;
    easing: string;
}
export interface CSSSKeyframeNode extends ASTNode {
    type: ASTNodeType.CSSS_KEYFRAME;
    percentage: number;
    declarations: CSSSDeclarationNode[];
}
export interface ShapeBehaviorNode extends ASTNode {
    type: ASTNodeType.SHAPE_BEHAVIOR;
    name: string;
    physics: ShapePhysicsNode[];
    events: ShapeEventNode[];
}
export interface ShapePhysicsNode extends ASTNode {
    type: ASTNodeType.SHAPE_PHYSICS;
    forces: ShapeForceNode[];
    constraints: ShapeConstraintNode[];
    matterState: MatterStateNode;
}
export interface ShapeForceNode extends ASTNode {
    type: ASTNodeType.SHAPE_FORCE;
    forceType: 'gravity' | 'elastic' | 'viscous' | 'electromagnetic';
    magnitude: ExpressionNode;
    direction: SphericalCoordinateNode;
}
export interface ShapeConstraintNode extends ASTNode {
    type: ASTNodeType.SHAPE_CONSTRAINT;
    constraintType: 'spherical_surface' | 'radial_range' | 'angular_cone';
    parameters: Record<string, ExpressionNode>;
}
export interface ShapeEventNode extends ASTNode {
    type: ASTNodeType.SHAPE_EVENT;
    trigger: string;
    response: ExpressionNode;
}
export interface StyleBotNode extends ASTNode {
    type: ASTNodeType.STYB_BOT;
    name: string;
    agents: StyleBotAgentNode[];
    parallel: boolean;
}
export interface StyleBotAgentNode extends ASTNode {
    type: ASTNodeType.STYB_AGENT;
    name: string;
    render: StyleBotRenderNode[];
    optimize: ExpressionNode[];
}
export interface StyleBotRenderNode extends ASTNode {
    type: ASTNodeType.STYB_RENDER;
    target: string;
    quality: ExpressionNode;
    priority: number;
}
export interface ExpressionNode extends ASTNode {
    type: ASTNodeType.EXPRESSION;
    expression: BinaryOperationNode | UnaryOperationNode | CallExpressionNode | MemberExpressionNode | LiteralNode | IdentifierNode;
}
export type ExpressionNodeType = BinaryOperationNode | UnaryOperationNode | CallExpressionNode | MemberExpressionNode | LiteralNode | IdentifierNode;
export interface BinaryOperationNode extends ASTNode {
    type: ASTNodeType.BINARY_OPERATION;
    operator: string;
    left: ExpressionNode;
    right: ExpressionNode;
}
export interface UnaryOperationNode extends ASTNode {
    type: ASTNodeType.UNARY_OPERATION;
    operator: string;
    operand: ExpressionNode;
}
export interface CallExpressionNode extends ASTNode {
    type: ASTNodeType.CALL_EXPRESSION;
    callee: ExpressionNode;
    arguments: ExpressionNode[];
}
export interface MemberExpressionNode extends ASTNode {
    type: ASTNodeType.MEMBER_EXPRESSION;
    object: ExpressionNode;
    property: string;
}
export interface LiteralNode extends ASTNode {
    type: ASTNodeType.LITERAL;
    value: string | number | boolean;
    raw: string;
}
export interface IdentifierNode extends ASTNode {
    type: ASTNodeType.IDENTIFIER;
    name: string;
}
export interface SphericalCoordinateNode extends ASTNode {
    type: ASTNodeType.SPHERICAL_COORDINATE;
    r: ExpressionNode;
    theta: ExpressionNode;
    phi: ExpressionNode;
}
export interface SolidAngleNode extends ASTNode {
    type: ASTNodeType.SOLID_ANGLE;
    omega: ExpressionNode;
    theta_min: ExpressionNode;
    theta_max: ExpressionNode;
    phi_min: ExpressionNode;
    phi_max: ExpressionNode;
}
export interface MatterStateNode extends ASTNode {
    type: ASTNodeType.MATTER_STATE;
    state: 'solid' | 'liquid' | 'gas' | 'plasma';
    properties: Record<string, ExpressionNode>;
}
export interface IfStatementNode extends ASTNode {
    type: ASTNodeType.IF_STATEMENT;
    condition: ExpressionNode;
    consequent: ASTNode;
    alternate?: ASTNode;
}
export interface WhileStatementNode extends ASTNode {
    type: ASTNodeType.WHILE_STATEMENT;
    condition: ExpressionNode;
    body: ASTNode;
}
export interface ForStatementNode extends ASTNode {
    type: ASTNodeType.FOR_STATEMENT;
    init?: ASTNode;
    condition?: ExpressionNode;
    update?: ASTNode;
    body: ASTNode;
}
export interface ReturnStatementNode extends ASTNode {
    type: ASTNodeType.RETURN_STATEMENT;
    argument?: ExpressionNode;
}
export declare class HSMLMultiParser {
    private tokens;
    private current;
    private language;
    constructor(tokens: Token[], language: 'hsml' | 'csss' | 'shape' | 'styb');
    parse(): ProgramNode;
    private parseTopLevel;
    private parseHSMLElement;
    private parseHSMLAttribute;
    private parseCSSSRule;
    private parseCSSSSelector;
    private parseCSSSDeclaration;
    private parseShapeBehavior;
    private parseShapePhysics;
    private parseShapeForce;
    private parseShapeConstraint;
    private parseShapeEvent;
    private parseStyleBot;
    private parseStyleBotAgent;
    private parseStyleBotRender;
    private parseExpression;
    private parseEquality;
    private parseComparison;
    private parseTerm;
    private parseFactor;
    private parseUnary;
    private parseCall;
    private finishCall;
    private parsePrimary;
    private parseSphericalPrimary;
    private parseSphericalCoordinate;
    private parseSolidAngle;
    private parseMatterState;
    private match;
    private check;
    private advance;
    private isAtEnd;
    private peek;
    private previous;
    private consume;
    private isMatterState;
    private synchronize;
    private parseDeclaration;
    private parseBlock;
    private parseIfStatement;
    private parseWhileStatement;
    private parseForStatement;
    private parseFunction;
    private parseReturnStatement;
}
//# sourceMappingURL=hsml-parser.d.ts.map