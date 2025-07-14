/**
 * HSML-SDT Development Server
 * ===========================
 * 
 * Integrated development server with hot reloading, physics simulation,
 * and zero-division safety monitoring
 */

import * as http from 'http';
import * as fs from 'fs';
import * as path from 'path';
import * as url from 'url';
import { WebSocketServer } from 'ws';

interface ServerConfig {
    port: number;
    host: string;
    open: boolean;
    watch: boolean;
    physics: boolean;
    collisions: boolean;
    root: string;
}

interface FileChange {
    type: 'added' | 'changed' | 'removed';
    path: string;
    timestamp: number;
}

interface PhysicsState {
    elements: number;
    collisions: number;
    frameRate: number;
    safetyViolations: number;
}

class HSMLSDTDevServer {
    private config: ServerConfig;
    private server: http.Server;
    private wss: WebSocketServer;
    private watcher: any = null;
    private clients: Set<any> = new Set();
    private physicsState: PhysicsState;

    constructor(config: ServerConfig) {
        this.config = config;
        this.server = http.createServer(this.handleRequest.bind(this));
        this.wss = new WebSocketServer({ server: this.server });
        this.physicsState = {
            elements: 0,
            collisions: 0,
            frameRate: 60,
            safetyViolations: 0
        };
        this.setupWebSocket();
    }

    public async start(): Promise<void> {
        try {
            console.log(`🌐 Starting HSML-SDT development server...`);
            console.log(`📍 Server: http://${this.config.host}:${this.config.port}`);
            console.log(`📁 Root: ${this.config.root}`);
            console.log(`⚡ Physics: ${this.config.physics ? '✅ Enabled' : '❌ Disabled'}`);
            console.log(`💥 Collisions: ${this.config.collisions ? '✅ Enabled' : '❌ Disabled'}`);
            
            // Start file watching if enabled
            if (this.config.watch) {
                await this.startFileWatching();
            }

            // Start physics monitoring if enabled
            if (this.config.physics) {
                this.startPhysicsMonitoring();
            }

            // Start the server
            this.server.listen(this.config.port, this.config.host, () => {
                console.log(`✅ HSML-SDT development server running`);
                console.log(`🔗 Open: http://${this.config.host}:${this.config.port}`);
                
                if (this.config.open) {
                    this.openBrowser();
                }
            });

        } catch (error) {
            console.error('❌ Failed to start HSML-SDT development server:', error);
            throw error;
        }
    }

    private setupWebSocket(): void {
        this.wss.on('connection', (ws) => {
            console.log('🔌 Client connected to HSML-SDT dev server');
            this.clients.add(ws);

            // Send initial configuration
            ws.send(JSON.stringify({
                type: 'config',
                data: {
                    physics: this.config.physics,
                    collisions: this.config.collisions,
                    features: ['21d-physics', 'zero-division-safety', 'hot-reload']
                }
            }));

            // Send current physics state
            if (this.config.physics) {
                ws.send(JSON.stringify({
                    type: 'physics-state',
                    data: this.physicsState
                }));
            }

            ws.on('close', () => {
                console.log('🔌 Client disconnected from HSML-SDT dev server');
                this.clients.delete(ws);
            });

            // Handle client messages
            ws.on('message', (message) => {
                try {
                    const data = JSON.parse(message.toString());
                    this.handleClientMessage(ws, data);
                } catch (error) {
                    console.error('Error parsing client message:', error);
                }
            });
        });
    }

    private handleClientMessage(ws: any, data: any): void {
        switch (data.type) {
            case 'physics-update':
                // Update physics state from client
                this.updatePhysicsState(data.state);
                break;
            
            case 'safety-violation':
                // Log safety violations
                this.physicsState.safetyViolations++;
                console.warn(`⚠️  Safety violation detected: ${data.details}`);
                this.broadcastToClients({
                    type: 'safety-alert',
                    data: { violation: data.details, count: this.physicsState.safetyViolations }
                });
                break;
            
            case 'performance-data':
                // Update performance metrics
                this.physicsState.frameRate = data.frameRate || 60;
                break;
            
            case 'request-reload':
                // Manual reload request
                this.broadcastReload();
                break;
        }
    }

    private async startFileWatching(): Promise<void> {
        try {
            // Use dynamic import for chokidar since it might not be available
            const chokidar = await import('chokidar');
            
            const watchPaths = [
                path.join(this.config.root, '**/*.ts'),
                path.join(this.config.root, '**/*.js'),
                path.join(this.config.root, '**/*.hsml'),
                path.join(this.config.root, '**/*.csss'),
                path.join(this.config.root, '**/*.html')
            ];

            this.watcher = chokidar.watch(watchPaths, {
                ignored: [
                    '**/node_modules/**',
                    '**/dist/**',
                    '**/.git/**'
                ],
                persistent: true,
                ignoreInitial: true
            });

            this.watcher.on('change', (filePath: string) => {
                console.log(`📝 File changed: ${filePath}`);
                this.handleFileChange({
                    type: 'changed',
                    path: filePath,
                    timestamp: Date.now()
                });
            });

            this.watcher.on('add', (filePath: string) => {
                console.log(`➕ File added: ${filePath}`);
                this.handleFileChange({
                    type: 'added',
                    path: filePath,
                    timestamp: Date.now()
                });
            });

            this.watcher.on('unlink', (filePath: string) => {
                console.log(`➖ File removed: ${filePath}`);
                this.handleFileChange({
                    type: 'removed',
                    path: filePath,
                    timestamp: Date.now()
                });
            });

            console.log('👀 File watching started');

        } catch (error) {
            console.warn('⚠️  File watching not available (chokidar not installed)');
            console.log('💡 Install chokidar for hot reloading: npm install chokidar');
        }
    }

    private handleFileChange(change: FileChange): void {
        // Determine if recompilation is needed
        const needsRecompile = this.needsRecompilation(change.path);
        
        if (needsRecompile) {
            this.recompileProject(change.path);
        }

        // Broadcast change to clients
        this.broadcastToClients({
            type: 'file-change',
            data: change
        });

        // If it's a source file, trigger reload
        if (this.shouldTriggerReload(change.path)) {
            setTimeout(() => this.broadcastReload(), 500); // Small delay for compilation
        }
    }

    private needsRecompilation(filePath: string): boolean {
        const ext = path.extname(filePath);
        return ['.ts', '.hsml', '.csss'].includes(ext);
    }

    private shouldTriggerReload(filePath: string): boolean {
        const ext = path.extname(filePath);
        return ['.ts', '.js', '.html', '.hsml', '.csss'].includes(ext);
    }

    private async recompileProject(changedFile: string): Promise<void> {
        try {
            console.log(`🔨 Recompiling project (${path.basename(changedFile)} changed)...`);
            
            // Run TypeScript compilation
            const { execSync } = require('child_process');
            execSync('npx tsc --build', { stdio: 'pipe' });
            
            console.log('✅ Recompilation completed');
            
        } catch (error) {
            console.error('❌ Recompilation failed:', error);
            
            // Send compilation error to clients
            this.broadcastToClients({
                type: 'compilation-error',
                data: { error: error.toString(), file: changedFile }
            });
        }
    }

    private startPhysicsMonitoring(): void {
        console.log('🔬 Starting physics monitoring...');
        
        // Monitor physics state every second
        setInterval(() => {
            this.broadcastToClients({
                type: 'physics-state',
                data: this.physicsState
            });
        }, 1000);

        // Monitor for safety violations
        setInterval(() => {
            if (this.physicsState.safetyViolations > 0) {
                console.log(`⚠️  Safety violations detected: ${this.physicsState.safetyViolations}`);
            }
        }, 5000);
    }

    private updatePhysicsState(newState: Partial<PhysicsState>): void {
        this.physicsState = { ...this.physicsState, ...newState };
    }

    private broadcastToClients(message: any): void {
        const messageStr = JSON.stringify(message);
        
        this.clients.forEach(client => {
            if (client.readyState === 1) { // WebSocket.OPEN
                client.send(messageStr);
            }
        });
    }

    private broadcastReload(): void {
        console.log('🔄 Broadcasting reload to clients...');
        
        this.broadcastToClients({
            type: 'reload',
            data: { timestamp: Date.now() }
        });
    }

    private handleRequest(req: http.IncomingMessage, res: http.ServerResponse): void {
        const parsedUrl = url.parse(req.url || '/', true);
        let pathname = parsedUrl.pathname || '/';

        // Normalize pathname
        if (pathname === '/') {
            pathname = '/index.html';
        }

        // Security: prevent directory traversal
        pathname = path.normalize(pathname);
        if (pathname.includes('..')) {
            res.writeHead(403);
            res.end('Forbidden');
            return;
        }

        // Try to serve the file
        this.serveFile(pathname, res);
    }

    private serveFile(pathname: string, res: http.ServerResponse): void {
        // Try multiple locations
        const possiblePaths = [
            path.join(this.config.root, pathname),
            path.join(this.config.root, 'demos', pathname),
            path.join(this.config.root, 'dist', pathname),
            path.join(this.config.root, 'public', pathname)
        ];

        let filePath: string | null = null;
        
        for (const possiblePath of possiblePaths) {
            if (fs.existsSync(possiblePath) && fs.statSync(possiblePath).isFile()) {
                filePath = possiblePath;
                break;
            }
        }

        // If no file found, try serving index.html from demos
        if (!filePath && pathname === '/index.html') {
            const demoIndex = path.join(this.config.root, 'demos', 'phase6-zero-division-fix-demo.html');
            if (fs.existsSync(demoIndex)) {
                filePath = demoIndex;
            }
        }

        if (!filePath) {
            this.serve404(res, pathname);
            return;
        }

        try {
            const ext = path.extname(filePath);
            const contentType = this.getContentType(ext);
            
            let content = fs.readFileSync(filePath, 'utf8');
            
            // Inject development features into HTML files
            if (ext === '.html') {
                content = this.injectDevFeatures(content);
            }

            res.writeHead(200, { 
                'Content-Type': contentType,
                'Cache-Control': 'no-cache'
            });
            res.end(content);
            
        } catch (error) {
            console.error('Error serving file:', error);
            res.writeHead(500);
            res.end('Internal Server Error');
        }
    }

    private injectDevFeatures(htmlContent: string): string {
        // Inject WebSocket client for hot reloading and physics monitoring
        const devScript = `
<script>
(function() {
    console.log('🔌 HSML-SDT Dev Server connecting...');
    
    const ws = new WebSocket('ws://${this.config.host}:${this.config.port}');
    let physicsEnabled = false;
    let collisionsEnabled = false;
    
    ws.onopen = function() {
        console.log('✅ Connected to HSML-SDT dev server');
    };
    
    ws.onmessage = function(event) {
        const message = JSON.parse(event.data);
        
        switch (message.type) {
            case 'config':
                physicsEnabled = message.data.physics;
                collisionsEnabled = message.data.collisions;
                console.log('⚙️  Dev server config:', message.data);
                break;
                
            case 'file-change':
                console.log('📝 File changed:', message.data.path);
                break;
                
            case 'reload':
                console.log('🔄 Reloading page...');
                window.location.reload();
                break;
                
            case 'compilation-error':
                console.error('❌ Compilation error:', message.data.error);
                showDevNotification('Compilation Error: ' + message.data.error, 'error');
                break;
                
            case 'physics-state':
                if (physicsEnabled) {
                    updatePhysicsDisplay(message.data);
                }
                break;
                
            case 'safety-alert':
                console.warn('⚠️  Safety violation:', message.data.violation);
                showDevNotification('Safety Violation: ' + message.data.violation, 'warning');
                break;
        }
    };
    
    ws.onclose = function() {
        console.log('🔌 Disconnected from dev server');
        showDevNotification('Dev server disconnected', 'warning');
    };
    
    // Function to show development notifications
    function showDevNotification(message, type = 'info') {
        const notification = document.createElement('div');
        notification.style.cssText = \`
            position: fixed;
            top: 20px;
            right: 20px;
            padding: 10px 20px;
            background: \${type === 'error' ? '#ff4444' : type === 'warning' ? '#ffaa44' : '#44ff44'};
            color: white;
            border-radius: 5px;
            font-family: monospace;
            z-index: 10000;
            max-width: 300px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.3);
        \`;
        notification.textContent = message;
        
        document.body.appendChild(notification);
        
        setTimeout(() => {
            notification.remove();
        }, 5000);
    }
    
    // Function to update physics display
    function updatePhysicsDisplay(state) {
        // Update or create physics display
        let display = document.getElementById('dev-physics-display');
        if (!display) {
            display = document.createElement('div');
            display.id = 'dev-physics-display';
            display.style.cssText = \`
                position: fixed;
                top: 20px;
                left: 20px;
                background: rgba(0,0,0,0.8);
                color: #00ff00;
                padding: 10px;
                border-radius: 5px;
                font-family: monospace;
                font-size: 12px;
                z-index: 9999;
                border: 1px solid #00ff00;
            \`;
            document.body.appendChild(display);
        }
        
        display.innerHTML = \`
            <div>📊 HSML-SDT Physics Monitor</div>
            <div>Elements: \${state.elements}</div>
            <div>Collisions: \${state.collisions}</div>
            <div>FPS: \${state.frameRate}</div>
            <div>Safety Violations: \${state.safetyViolations}</div>
        \`;
    }
    
    // Report physics updates to server
    function reportPhysicsUpdate(state) {
        if (ws.readyState === WebSocket.OPEN) {
            ws.send(JSON.stringify({
                type: 'physics-update',
                state: state
            }));
        }
    }
    
    // Report safety violations
    function reportSafetyViolation(details) {
        if (ws.readyState === WebSocket.OPEN) {
            ws.send(JSON.stringify({
                type: 'safety-violation',
                details: details
            }));
        }
    }
    
    // Make functions available globally
    window.devServer = {
        reportPhysicsUpdate,
        reportSafetyViolation
    };
    
})();
</script>`;

        // Inject before closing body tag, or at the end if no body tag
        if (htmlContent.includes('</body>')) {
            return htmlContent.replace('</body>', devScript + '</body>');
        } else {
            return htmlContent + devScript;
        }
    }

    private getContentType(ext: string): string {
        const mimeTypes: { [key: string]: string } = {
            '.html': 'text/html',
            '.js': 'text/javascript',
            '.css': 'text/css',
            '.json': 'application/json',
            '.png': 'image/png',
            '.jpg': 'image/jpeg',
            '.gif': 'image/gif',
            '.svg': 'image/svg+xml',
            '.ico': 'image/x-icon',
            '.hsml': 'text/plain',
            '.csss': 'text/plain'
        };
        
        return mimeTypes[ext] || 'text/plain';
    }

    private serve404(res: http.ServerResponse, pathname: string): void {
        const html404 = `
<!DOCTYPE html>
<html>
<head>
    <title>404 - Not Found | HSML-SDT Dev Server</title>
    <style>
        body { 
            font-family: monospace; 
            background: #0a0a0a; 
            color: #00ff00; 
            padding: 50px;
            text-align: center;
        }
        .error { color: #ff4444; }
        .suggestion { margin: 20px 0; padding: 20px; background: #001100; border-radius: 5px; }
    </style>
</head>
<body>
    <h1 class="error">404 - Not Found</h1>
    <p>File not found: <code>${pathname}</code></p>
    
    <div class="suggestion">
        <h3>🎮 Available Demos:</h3>
        <ul style="list-style: none; padding: 0;">
            <li><a href="/phase6-zero-division-fix-demo.html" style="color: #00ff00;">Phase 6 Zero-Division Fix Demo</a></li>
            <li><a href="/authentic-21d-demo.html" style="color: #00ff00;">Authentic 21D Framework Demo</a></li>
            <li><a href="/revolutionary-sdt-demo.html" style="color: #00ff00;">Revolutionary SDT Demo</a></li>
        </ul>
    </div>
    
    <div class="suggestion">
        <h3>🔧 Dev Server Features:</h3>
        <ul style="list-style: none; padding: 0;">
            <li>✅ Hot Reloading</li>
            <li>✅ Physics Monitoring</li>
            <li>✅ Zero-Division Safety Alerts</li>
            <li>✅ Compilation Error Display</li>
        </ul>
    </div>
</body>
</html>`;
        
        res.writeHead(404, { 'Content-Type': 'text/html' });
        res.end(html404);
    }

    private openBrowser(): void {
        const url = `http://${this.config.host}:${this.config.port}`;
        
        try {
            const { execSync } = require('child_process');
            const platform = process.platform;
            
            if (platform === 'darwin') {
                execSync(`open ${url}`);
            } else if (platform === 'win32') {
                execSync(`start ${url}`);
            } else {
                execSync(`xdg-open ${url}`);
            }
            
        } catch (error) {
            console.log(`🌐 Please open: ${url}`);
        }
    }

    public stop(): void {
        if (this.watcher) {
            this.watcher.close();
        }
        
        this.server.close();
        this.wss.close();
        
        console.log('🛑 HSML-SDT development server stopped');
    }
}

// Command line interface
if (require.main === module) {
    const args = process.argv.slice(2);
    
    const config: ServerConfig = {
        port: 3000,
        host: 'localhost',
        open: false,
        watch: true,
        physics: true,
        collisions: true,
        root: process.cwd()
    };
    
    // Parse command line arguments
    for (let i = 0; i < args.length; i += 2) {
        const arg = args[i];
        const value = args[i + 1];
        
        switch (arg) {
            case '--port':
                config.port = parseInt(value);
                break;
            case '--host':
                config.host = value;
                break;
            case '--open':
                config.open = true;
                i--; // No value for this flag
                break;
            case '--watch':
                config.watch = value !== 'false';
                break;
            case '--physics':
                config.physics = value !== 'false';
                break;
            case '--collisions':
                config.collisions = value !== 'false';
                break;
        }
    }
    
    const server = new HSMLSDTDevServer(config);
    server.start().catch(error => {
        console.error('Failed to start development server:', error);
        process.exit(1);
    });
    
    // Graceful shutdown
    process.on('SIGINT', () => {
        console.log('\n🛑 Shutting down HSML-SDT development server...');
        server.stop();
        process.exit(0);
    });
}

export { HSMLSDTDevServer };