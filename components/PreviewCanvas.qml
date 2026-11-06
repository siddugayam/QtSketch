import QtQuick
import QtQuick.Controls
import QtQuick.Layouts

// Simple preview canvas for Settings dialog
// Shows the molecule with draft settings
Canvas {
    id: root
    
    property var primitives: null
    property real draftBondLength: 1.5
    property real draftLineWidth: 0.6
    property string draftFontFamily: "Arial"
    property int draftFontSize: 10
    property bool draftMonochrome: true
    
    property real scale: 1.0
    property real offsetX: 0
    property real offsetY: 0
    
    antialiasing: true
    
    onPaint: {
        const ctx = getContext("2d")
        ctx.clearRect(0, 0, width, height)
        
        if (!primitives || !primitives.atoms || !primitives.bonds) return
        
        // Center the molecule
        const margin = 20
        const centerX = width / 2
        const centerY = height / 2
        
        // Draw bonds
        for (let i = 0; i < primitives.bonds.length; i++) {
            const b = primitives.bonds[i]
            const a1 = primitives.atomsById[b.begin]
            const a2 = primitives.atomsById[b.end]
            if (!a1 || !a2) continue
            
            // Scale coordinates
            const x1 = centerX + a1.x * scale * 40
            const y1 = centerY + a1.y * scale * 40
            const x2 = centerX + a2.x * scale * 40
            const y2 = centerY + a2.y * scale * 40
            
            // Set line width based on draft settings
            ctx.lineWidth = draftLineWidth * scale
            ctx.strokeStyle = draftMonochrome ? "#000000" : "#000000"
            
            ctx.beginPath()
            ctx.moveTo(x1, y1)
            ctx.lineTo(x2, y2)
            ctx.stroke()
        }
        
        // Draw atoms
        for (let i = 0; i < primitives.atoms.length; i++) {
            const a = primitives.atoms[i]
            const x = centerX + a.x * scale * 40
            const y = centerY + a.y * scale * 40
            
            // Set font based on draft settings
            ctx.font = (a.bold ? "bold " : "") + (a.italic ? "italic " : "") + draftFontSize + "px " + draftFontFamily
            ctx.fillStyle = draftMonochrome ? "#000000" : "#000000"
            ctx.textAlign = "center"
            ctx.textBaseline = "middle"
            
            // Draw atom label
            let label = a.label
            if (!label && a.element) label = a.element
            if (label) {
                ctx.fillText(label, x, y)
            }
        }
    }
    
    onWidthChanged: requestPaint()
    onHeightChanged: requestPaint()
    onPrimitivesChanged: requestPaint()
}
