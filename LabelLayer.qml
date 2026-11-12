import QtQuick
import QtQuick.Controls
Item {
    id: root
    anchors.fill: parent

    property real scale: 1.0
    property real offsetX: 0
    property real offsetY: 0
    property real bondLength: 40
    property Item canvas
    property int renderVersion: 0

    // Populated by labelCanvas.onPaint each repaint -- one entry per drawn
    // check-warning badge, since the badge itself is raster-painted rather
    // than a discrete hoverable Item. checkHoverArea below hit-tests against
    // this list to show the underlying issue text as a tooltip.
    property var checkWarningHotspots: []

    // a.checkWarning carries Indigo's raw check-type code (see IndigoService.cpp's
    // atomChecks set) -- readable text for the ones we know; anything else falls
    // back to the code itself with underscores turned to spaces rather than
    // showing nothing.
    readonly property var _checkTypeLabels: ({
        valence: "Valence error",
        radical: "Unusual radical state",
        pseudoatom: "Pseudoatom",
        stereo: "Ambiguous stereochemistry",
        ambiguous_h: "Ambiguous hydrogen count",
        "3d_coord": "Non-planar (3D) coordinates",
        overlap_atom: "Overlapping atoms",
        overlap_bond: "Overlapping bonds"
    })
    function checkTypeLabel(code) {
        return root._checkTypeLabels[code] || code.replace(/_/g, " ")
    }

    // ctx.measureText is re-run for every label component on every repaint; cache by
    // (font, text) since results don't change unless the font string does (scale/style).
    property var _measureCache: ({})
    function measureCached(ctx, text) {
        const key = ctx.font + "|" + text
        var v = _measureCache[key]
        if (v === undefined) {
            if (Object.keys(_measureCache).length > 5000) _measureCache = ({})
            v = ctx.measureText(text).width
            _measureCache[key] = v
        }
        return v
    }

    Canvas {
        id: labelCanvas
        anchors.fill: parent
        antialiasing: true

        onWidthChanged: requestPaint()
        onHeightChanged: requestPaint()

        function drawRoundedRect(ctx, x, y, w, h, r) {
            ctx.moveTo(x + r, y)
            ctx.lineTo(x + w - r, y)
            ctx.arcTo(x + w, y, x + w, y + r, r)
            ctx.lineTo(x + w, y + h - r)
            ctx.arcTo(x + w, y + h, x + w - r, y + h, r)
            ctx.lineTo(x + r, y + h)
            ctx.arcTo(x, y + h, x, y + h - r, r)
            ctx.lineTo(x, y + r)
            ctx.arcTo(x, y, x + r, y, r)
        }

        onPaint: {
            const ctx = getContext("2d")
            ctx.clearRect(0, 0, width, height)

            const newHotspots = []
            root.checkWarningHotspots = newHotspots

            if (!canvas) return
            if (!canvas.sketch.primitives || !canvas.sketch.primitives.atoms) return

            const fontSize = Theme.clampFontSize(Theme.baseFontSize, root.scale)
            const superSize = Math.max(7, 10 * root.scale)
            const subSize   = Math.max(7, 10 * root.scale)

            // Cheap pre-pass: approximate canvas position + "occupied radius" for
            // every atom, used only to steer check-warning badges away from
            // overlapping ANY atom's own label pill (including the flagged atom's
            // own pill, and any neighbour's). A generous constant multiple of
            // pillH stands in for the real per-atom pill width -- exact text
            // measurement isn't needed for a collision-avoidance estimate, and
            // redoing it here for every atom would double the paint-time cost.
            const approxPillH = fontSize * 1.1 + 8 * root.scale
            const atomOccupancy = []
            for (let ai = 0; ai < canvas.sketch.primitives.atoms.length; ai++) {
                const aa = canvas.sketch.primitives.atoms[ai]
                if (aa.isSgroup) { atomOccupancy.push(null); continue }
                const hasVisibleLabel = aa.isSgroup || aa.element !== "C" || aa.charge !== 0 ||
                        (aa.cipLabel && aa.cipLabel !== "") || canvas.showExplicitH ||
                        (aa.isotope && aa.isotope > 0) || (aa.radical && aa.radical > 0) ||
                        (aa.explicitValence !== undefined && aa.explicitValence >= 0) ||
                        (aa.attachmentPoints && aa.attachmentPoints > 0) ||
                        (aa.implicitHCount || 0) > 0
                if (!hasVisibleLabel) { atomOccupancy.push(null); continue }
                const ap = canvas.chemToCanvas(aa.x, aa.y)
                atomOccupancy.push({ x: ap.x, y: ap.y, r: approxPillH * 0.9 })
            }

            // Chiral-center CIP labels (R/S) are drawn beside their own wedge/hash
            // bond rather than centered on the atom -- map atomId -> unit direction
            // of that bond (the stereocenter is always the narrow/begin end).
            const chiralBondDir = {}
            if (canvas.sketch.primitives.bonds && canvas.sketch.primitives.atomsById) {
                const bondsList = canvas.sketch.primitives.bonds
                for (let bi = 0; bi < bondsList.length; bi++) {
                    const bb = bondsList[bi]
                    if (bb.stereo !== 1 && bb.stereo !== 6) continue
                    if (chiralBondDir[bb.begin] !== undefined) continue
                    const beginA = canvas.sketch.primitives.atomsById[bb.begin.toString()]
                    const endA = canvas.sketch.primitives.atomsById[bb.end.toString()]
                    if (!beginA || !endA) continue
                    const bp1 = canvas.chemToCanvas(beginA.x, beginA.y)
                    const bp2 = canvas.chemToCanvas(endA.x, endA.y)
                    const bdx = bp2.x - bp1.x, bdy = bp2.y - bp1.y
                    const blen = Math.sqrt(bdx * bdx + bdy * bdy)
                    if (blen > 0) chiralBondDir[bb.begin] = { x: bdx / blen, y: bdy / blen }
                }
            }

            // Given a flagged atom's own canvas position and pill half-extents,
            // pick the closest-to-original badge center (searching 8 compass
            // directions, preferring the original upper-left placement) that
            // does not overlap any atom's own label pill. Falls back to the
            // original upper-left placement if every direction collides.
            function pickBadgeCenter(px, py, pillHalfW, pillHalfH, warnR, selfIdx) {
                const baseDist = Math.max(pillHalfW, pillHalfH) + warnR + 2 * root.scale
                const angles = [225, 180, 270, 135, 315, 90, 0, 45] // degrees; 225 = up-left (original)
                let fallback = null
                for (let k = 0; k < angles.length; k++) {
                    const rad = angles[k] * Math.PI / 180
                    const cx = px + baseDist * Math.cos(rad)
                    const cy = py + baseDist * Math.sin(rad)
                    if (k === 0) fallback = { x: cx, y: cy }
                    let collides = false
                    for (let oi = 0; oi < atomOccupancy.length; oi++) {
                        if (oi === selfIdx) continue
                        const occ = atomOccupancy[oi]
                        if (!occ) continue
                        const ddx = cx - occ.x, ddy = cy - occ.y
                        const minDist = occ.r + warnR
                        if (ddx * ddx + ddy * ddy < minDist * minDist) { collides = true; break }
                    }
                    if (!collides) return { x: cx, y: cy }
                }
                return fallback
            }

            for (let i = 0; i < canvas.sketch.primitives.atoms.length; i++) {
                const a = canvas.sketch.primitives.atoms[i]

                // ── Contracted sgroup node (FG / S&S abbreviation pill) ──────
                if (a.isSgroup) {
                    const p = canvas.chemToCanvas(a.x, a.y)
                    const isSelected = canvas.sketch.selection.atom_ids && canvas.sketch.selection.atom_ids.indexOf(a.id) >= 0
                    ctx.font = "bold " + fontSize + "px " + Theme.fontFamilyCss
                    const tw = root.measureCached(ctx, a.label)
                    const padX = 5 * root.scale
                    const padY = 3 * root.scale
                    const pillW = tw + padX * 2
                    const pillH = fontSize + padY * 2
                    const rx = pillH / 2
                    // Background pill
                    ctx.fillStyle = isSelected ? Theme.selected : Theme.background
                    ctx.strokeStyle = isSelected ? Theme.accent : Theme.textPrimary
                    ctx.lineWidth = Theme.clampStrokeWidth(1.2, root.scale)
                    ctx.beginPath()
                    drawRoundedRect(ctx, p.x - pillW / 2, p.y - pillH / 2, pillW, pillH, rx)
                    ctx.fill()
                    ctx.stroke()
                    // Label text
                    ctx.fillStyle = isSelected ? Theme.accent : Theme.textPrimary
                    ctx.textAlign = "center"
                    ctx.textBaseline = "middle"
                    ctx.fillText(a.label, p.x, p.y)
                    continue
                }
                // ────────────────────────────────────────────────────────────

                const hasIsotope  = a.isotope && a.isotope > 0
                const hasRadical  = a.radical && a.radical > 0
                const hasValence  = a.explicitValence !== undefined && a.explicitValence >= 0
                const hasH        = (a.implicitHCount || 0) > 0
                const hasAttachment = a.attachmentPoints && a.attachmentPoints > 0

                // Skip unlabelled carbons
                if (a.element === "C" && a.charge === 0 && !a.cipLabel &&
                        !canvas.showExplicitH && !hasIsotope && !hasRadical && !hasValence && !hasAttachment) {
                    if (!hasH) continue
                    // Terminal carbons with showExplicitH off still show no label if renderLabel === "C"
                    if (a.label === "C") continue
                }

                const p = canvas.chemToCanvas(a.x, a.y)
                const isSelected = canvas.sketch.selection.atom_ids && canvas.sketch.selection.atom_ids.indexOf(a.id) >= 0

                const baseElement  = a.element
                const implicitH    = a.implicitHCount || 0
                const charge       = a.charge || 0
                const cipLabel     = a.cipLabel || ""
                const aam          = a.aam || 0
                const checkWarning = a.checkWarning || ""
                const hOnLeft      = a.hOnLeft || false

                // Chiral centers on a plain carbon show only the CIP letter (R/S/E/Z),
                // not the "C" element symbol -- matches standard skeletal-formula
                // convention (bare vertex + stereo descriptor, no atom label).
                const isChiralCarbon = baseElement === "C" && charge === 0 && !hasIsotope &&
                        !hasRadical && !hasValence && !hasAttachment && cipLabel !== ""

                // When nothing else needs a pill at this vertex, draw the CIP letter
                // as a small floating label beside the atom's own wedge/hash bond
                // instead of centered on the vertex -- keeps the bond/vertex clear.
                if (isChiralCarbon && implicitH === 0 && !canvas.showExplicitH && chiralBondDir[a.id]) {
                    const dir = chiralBondDir[a.id]
                    const perpX = -dir.y, perpY = dir.x
                    const cipFontSize = Math.max(9, 11 * root.scale)
                    const cipOffset = 10 * root.scale
                    ctx.font = "bold " + cipFontSize + "px " + Theme.fontFamilyCss
                    ctx.textAlign = "center"
                    ctx.textBaseline = "middle"
                    ctx.fillStyle = Theme.accent
                    ctx.fillText(cipLabel, p.x + perpX * cipOffset, p.y + perpY * cipOffset)
                    continue
                }

                // ── Measure components ──────────────────────────────────────
                ctx.font = "bold " + fontSize + "px " + Theme.fontFamilyCss
                const elemW = isChiralCarbon ? 0 : root.measureCached(ctx, baseElement)

                let hText = "", hSubText = ""
                if (implicitH > 0) {
                    hText = "H"
                    if (implicitH > 1) hSubText = implicitH.toString()
                }
                ctx.font = "bold " + fontSize + "px " + Theme.fontFamilyCss
                const hW = hText ? root.measureCached(ctx, hText) : 0
                ctx.font = "bold " + subSize + "px " + Theme.fontFamilyCss
                const hSubW = hSubText ? root.measureCached(ctx, hSubText) : 0

                let chargeText = ""
                if (charge !== 0) {
                    if      (charge ===  1) chargeText = "+"
                    else if (charge === -1) chargeText = "–"
                    else if (charge >   1)  chargeText = charge + "+"
                    else                    chargeText = Math.abs(charge) + "–"
                }
                ctx.font = "bold " + superSize + "px " + Theme.fontFamilyCss
                const chargeW = chargeText ? root.measureCached(ctx, chargeText) : 0

                // Attachment point badge ("*1", "*2" etc.)
                let apText = ""
                if (hasAttachment) {
                    apText = "*" + a.attachmentPoints
                }
                ctx.font = "bold " + superSize + "px " + Theme.fontFamilyCss
                const apW = apText ? root.measureCached(ctx, apText) : 0

                // Isotope prefix ("13" for ¹³C)
                const isoText = hasIsotope ? a.isotope.toString() : ""
                ctx.font = "bold " + superSize + "px " + Theme.fontFamilyCss
                const isoW = isoText ? root.measureCached(ctx, isoText) : 0

                // Explicit valence suffix ("(4)" etc.)
                const valText = hasValence ? ("(" + a.explicitValence + ")") : ""
                ctx.font = "bold " + subSize + "px " + Theme.fontFamilyCss
                const valW = valText ? root.measureCached(ctx, valText) : 0

                let totalWidth = isoW + elemW + hW + hSubW + chargeW + valW + apW
                let cipW = 0
                const cipText = isChiralCarbon ? cipLabel : " " + cipLabel
                if (cipLabel) {
                    ctx.font = "bold " + subSize + "px " + Theme.fontFamilyCss
                    cipW = root.measureCached(ctx, cipText)
                    totalWidth += cipW
                }

                // ── Background pill ──────────────────────────────────────────
                const padX = 4 * root.scale
                const padY = 4 * root.scale
                const pillW = totalWidth + padX * 2
                const pillH = fontSize * 1.1 + padY * 2

                ctx.fillStyle = Theme.surface
                ctx.beginPath()
                ctx.ellipse(p.x, p.y, pillW / 2, pillH / 2, 0, 0, Math.PI * 2)
                ctx.fill()

                // a.color is always set, so check the CPK switch here too.
                const colorStr = Theme.cpkColorsEnabled ? (a.color || Theme.getElementColor(a.element)) : Theme.textPrimary
                const textColor = isSelected ? Theme.accent : colorStr

                // ── Draw label left-to-right ─────────────────────────────────
                let cx = p.x - totalWidth / 2

                // Isotope superscript before element
                if (isoText) {
                    ctx.font = "bold " + superSize + "px " + Theme.fontFamilyCss
                    ctx.fillStyle = textColor
                    ctx.textAlign = "left"
                    ctx.textBaseline = "middle"
                    ctx.fillText(isoText, cx, p.y - fontSize * Theme.superscriptOffset)
                    cx += isoW
                }

                if (hOnLeft && hText) {
                    // H on left side
                    ctx.font = "bold " + fontSize + "px " + Theme.fontFamilyCss
                    ctx.fillStyle = textColor
                    ctx.textAlign = "left"
                    ctx.textBaseline = "middle"
                    ctx.fillText(hText, cx, p.y)
                    cx += hW
                    if (hSubText) {
                        ctx.font = "bold " + subSize + "px " + Theme.fontFamilyCss
                        ctx.fillText(hSubText, cx, p.y + fontSize * Theme.subscriptOffset)
                        cx += hSubW
                    }
                    if (!isChiralCarbon) {
                        ctx.font = "bold " + fontSize + "px " + Theme.fontFamilyCss
                        ctx.fillText(baseElement, cx, p.y)
                    }
                    cx += elemW
                } else {
                    // Element first, H on right
                    ctx.font = "bold " + fontSize + "px " + Theme.fontFamilyCss
                    ctx.fillStyle = textColor
                    ctx.textAlign = "left"
                    ctx.textBaseline = "middle"
                    if (!isChiralCarbon) {
                        ctx.fillText(baseElement, cx, p.y)
                    }
                    cx += elemW
                    if (hText) {
                        ctx.fillText(hText, cx, p.y)
                        cx += hW
                        if (hSubText) {
                            ctx.font = "bold " + subSize + "px " + Theme.fontFamilyCss
                            ctx.fillText(hSubText, cx, p.y + fontSize * Theme.subscriptOffset)
                            cx += hSubW
                        }
                    }
                }

                // Charge superscript
                if (chargeText) {
                    ctx.font = "bold " + superSize + "px " + Theme.fontFamilyCss
                    ctx.fillStyle = textColor
                    ctx.fillText(chargeText, cx, p.y - fontSize * Theme.superscriptOffset)
                    cx += chargeW
                }

                // Attachment point superscript
                if (apText) {
                    ctx.font = "bold " + superSize + "px " + Theme.fontFamilyCss
                    ctx.fillStyle = textColor
                    ctx.fillText(apText, cx, p.y - fontSize * Theme.superscriptOffset)
                    cx += apW
                }

                // Explicit valence subscript
                if (valText) {
                    ctx.font = "bold " + subSize + "px " + Theme.fontFamilyCss
                    ctx.fillStyle = isSelected ? Theme.accent : Theme.textSecondary
                    ctx.fillText(valText, cx, p.y + fontSize * Theme.subscriptOffset)
                    cx += valW
                }

                // CIP descriptor (R/S/E/Z)
                if (cipLabel) {
                    ctx.font = "bold " + subSize + "px " + Theme.fontFamilyCss
                    ctx.fillStyle = Theme.accent
                    ctx.fillText(cipText, cx, p.y)
                    cx += cipW
                }

                // ── Radical dots above pill ──────────────────────────────────
                if (hasRadical) {
                    const dotR  = Math.max(2, 2.5 * root.scale)
                    const dotY  = p.y - pillH / 2 - dotR - 1 * root.scale
                    ctx.fillStyle = textColor
                    if (a.radical === 1) {
                        // Monoradical: single dot
                        ctx.beginPath()
                        ctx.arc(p.x, dotY, dotR, 0, Math.PI * 2)
                        ctx.fill()
                    } else if (a.radical === 2) {
                        // Singlet diradical: two dots close together (paired)
                        const sep = dotR * 1.8
                        ctx.beginPath()
                        ctx.arc(p.x - sep / 2, dotY, dotR, 0, Math.PI * 2)
                        ctx.fill()
                        ctx.beginPath()
                        ctx.arc(p.x + sep / 2, dotY, dotR, 0, Math.PI * 2)
                        ctx.fill()
                        // Small line connecting them to signal pairing
                        ctx.strokeStyle = textColor
                        ctx.lineWidth = root.scale
                        ctx.beginPath()
                        ctx.moveTo(p.x - sep / 2 + dotR, dotY)
                        ctx.lineTo(p.x + sep / 2 - dotR, dotY)
                        ctx.stroke()
                    } else {
                        // Triplet diradical: two spaced dots
                        const sep = dotR * 4
                        ctx.beginPath()
                        ctx.arc(p.x - sep / 2, dotY, dotR, 0, Math.PI * 2)
                        ctx.fill()
                        ctx.beginPath()
                        ctx.arc(p.x + sep / 2, dotY, dotR, 0, Math.PI * 2)
                        ctx.fill()
                    }
                }

                // AAM (atom-atom mapping) number badge — small orange circle above pill
                if (aam > 0) {
                    const badgeR = Math.max(6, 8 * root.scale)
                    const badgeY = p.y - pillH / 2 - badgeR - 1 * root.scale
                    ctx.fillStyle = Theme.badgeAam
                    ctx.beginPath()
                    ctx.arc(p.x + pillW / 2 + badgeR, badgeY, badgeR, 0, Math.PI * 2)
                    ctx.fill()
                    ctx.fillStyle = Theme.badgeText
                    ctx.font = "bold " + Math.max(7, 9 * root.scale) + "px " + Theme.fontFamilyCss
                    ctx.textAlign = "center"
                    ctx.textBaseline = "middle"
                    ctx.fillText(aam.toString(), p.x + pillW / 2 + badgeR, badgeY)
                }

                // Check warning badge — small red circle with "!", steered away from
                // overlapping any atom's own label pill (see pickBadgeCenter above).
                if (checkWarning) {
                    const warnR = Math.max(5, 7 * root.scale)
                    const warnCenter = pickBadgeCenter(p.x, p.y, pillW / 2, pillH / 2, warnR, i)
                    const warnCx = warnCenter.x
                    const warnY = warnCenter.y
                    ctx.fillStyle = Theme.error
                    ctx.beginPath()
                    ctx.arc(warnCx, warnY, warnR, 0, Math.PI * 2)
                    ctx.fill()
                    ctx.fillStyle = Theme.badgeText
                    ctx.font = "bold " + Math.max(8, 10 * root.scale) + "px " + Theme.fontFamilyCss
                    ctx.textAlign = "center"
                    ctx.textBaseline = "middle"
                    ctx.fillText("!", warnCx, warnY)
                    // Recorded for the hover overlay below -- the badge is raster-
                    // painted, not a discrete Item, so hit-testing for the tooltip
                    // needs its own copy of these coordinates.
                    newHotspots.push({ x: warnCx, y: warnY, r: warnR + 2 * root.scale, text: root.checkTypeLabel(checkWarning) })
                }
            }

            // ── Text annotations ─────────────────────────────────────────────
            if (canvas.sketch.primitives.texts) {
                const txFont = Theme.clampFontSize(Theme.baseFontSize, root.scale)
                ctx.textAlign = "left"
                ctx.textBaseline = "top"
                ctx.fillStyle = Theme.textPrimary
                for (let ti = 0; ti < canvas.sketch.primitives.texts.length; ti++) {
                    const t = canvas.sketch.primitives.texts[ti]
                    const tp = canvas.chemToCanvas(t.x, t.y)
                    let fontStr = ""
                    if (t.italic) fontStr += "italic "
                    if (t.bold) fontStr += "bold "
                    fontStr += txFont + "px " + Theme.fontFamilyCss
                    ctx.font = fontStr

                    const txLines = String(t.content).split("\n")
                    for (let li = 0; li < txLines.length; li++) {
                        ctx.fillText(txLines[li], tp.x, tp.y + li * txFont * 1.25)
                    }
                }
            }

            // ── Bracket-around-group annotations ──────────────────────────────
            if (canvas.sketch.primitives.brackets) {
                ctx.strokeStyle = Theme.textPrimary
                ctx.lineWidth = 1.5
                for (let bi = 0; bi < canvas.sketch.primitives.brackets.length; bi++) {
                    const br = canvas.sketch.primitives.brackets[bi]
                    const topLeft = canvas.chemToCanvas(br.minX, br.minY)
                    const bottomRight = canvas.chemToCanvas(br.maxX, br.maxY)
                    const left = Math.min(topLeft.x, bottomRight.x)
                    const right = Math.max(topLeft.x, bottomRight.x)
                    const top = Math.min(topLeft.y, bottomRight.y)
                    const bottom = Math.max(topLeft.y, bottomRight.y)
                    const tick = 6 * root.scale

                    ctx.beginPath()
                    ctx.moveTo(left + tick, top)
                    ctx.lineTo(left, top)
                    ctx.lineTo(left, bottom)
                    ctx.lineTo(left + tick, bottom)
                    ctx.stroke()

                    ctx.beginPath()
                    ctx.moveTo(right - tick, top)
                    ctx.lineTo(right, top)
                    ctx.lineTo(right, bottom)
                    ctx.lineTo(right - tick, bottom)
                    ctx.stroke()
                }
            }
            // Bond CIP labels
            if (canvas.sketch.primitives.bonds && canvas.sketch.primitives.atomsById) {
                const bondsList = canvas.sketch.primitives.bonds
                for (let bi = 0; bi < bondsList.length; bi++) {
                    const bb = bondsList[bi]
                    if (!bb.cipLabel) continue

                    const beginA = canvas.sketch.primitives.atomsById[bb.begin.toString()]
                    const endA = canvas.sketch.primitives.atomsById[bb.end.toString()]
                    if (!beginA || !endA) continue

                    const p1 = canvas.chemToCanvas(beginA.x, beginA.y)
                    const p2 = canvas.chemToCanvas(endA.x, endA.y)

                    const midX = (p1.x + p2.x) / 2
                    const midY = (p1.y + p2.y) / 2

                    const bdx = p2.x - p1.x
                    const bdy = p2.y - p1.y
                    const blen = Math.sqrt(bdx * bdx + bdy * bdy)
                    if (blen === 0) continue

                    const dirX = bdx / blen
                    const dirY = bdy / blen
                    const perpX = -dirY
                    const perpY = dirX

                    const cipFontSize = Math.max(9, 11 * root.scale)
                    const cipOffset = 10 * root.scale

                    ctx.font = "bold " + cipFontSize + "px " + Theme.fontFamilyCss
                    ctx.textAlign = "center"
                    ctx.textBaseline = "middle"
                    ctx.fillStyle = Theme.accent
                    ctx.fillText(bb.cipLabel, midX + perpX * cipOffset, midY + perpY * cipOffset)
                }
            }

            // Images
            if (canvas.sketch.primitives.images) {
                for (let ii = 0; ii < canvas.sketch.primitives.images.length; ii++) {
                    const imgNode = canvas.sketch.primitives.images[ii]
                    let tp = canvas.chemToCanvas(imgNode.x, imgNode.y)
                    let pixelW = imgNode.w * canvas.chemScale
                    let pixelH = imgNode.h * canvas.chemScale

                    if (canvas.selectedImageId === imgNode.id) {
                        const ma = canvas.mouseArea
                        if (ma && ma.resizingSelection) {
                            const ax = ma.resizeAnchorCanvasX
                            const ay = ma.resizeAnchorCanvasY
                            const x1 = ax + (tp.x - ax) * ma.currentResizeFactor
                            const y1 = ay + (tp.y - ay) * ma.currentResizeFactor
                            const x2 = ax + (tp.x + pixelW - ax) * ma.currentResizeFactor
                            const y2 = ay + (tp.y + pixelH - ay) * ma.currentResizeFactor
                            tp = Qt.point(Math.min(x1, x2), Math.min(y1, y2))
                            pixelW = Math.abs(x2 - x1)
                            pixelH = Math.abs(y2 - y1)
                        } else if (ma && ma.movingImage) {
                            const dx = ma.mouseX - ma.pressX
                            const dy = ma.mouseY - ma.pressY
                            tp = Qt.point(tp.x + dx, tp.y + dy)
                        }
                    }

                    const imgItem = imageRepeater.itemAt(ii)
                    if (imgItem && imgItem.status === Image.Ready) {
                        ctx.drawImage(imgItem, tp.x, tp.y, pixelW, pixelH)
                    }
                }
            }
        }

        Connections {
            target: Theme
            function onCpkColorsEnabledChanged() { labelCanvas.requestPaint() }
        }

        Connections {
            target: root
            function onRenderVersionChanged() { labelCanvas.requestPaint() }
            function onScaleChanged()         { labelCanvas.requestPaint() }
            function onOffsetXChanged()        { if (!canvas || !canvas._panningActive) labelCanvas.requestPaint() }
            function onOffsetYChanged()        { if (!canvas || !canvas._panningActive) labelCanvas.requestPaint() }
        }
    }

    Repeater {
        id: imageRepeater
        model: canvas && canvas.sketch && canvas.sketch.primitives ? (canvas.sketch.primitives.images || []) : []
        Image {
            visible: false
            source: modelData.bitmap
            onStatusChanged: {
                if (status === Image.Ready)
                    labelCanvas.requestPaint()
                else if (status === Image.Error)
                    console.warn("Failed to load embedded image:", modelData.bitmap ? modelData.bitmap.substring(0, 40) : "")
            }
        }
    }

    // Hover-only overlay for check-warning badges: acceptedButtons is NoButton
    // so clicks/drags to atoms and bonds underneath are never intercepted --
    // hover tracking works independently of that in MouseArea.
    MouseArea {
        anchors.fill: parent
        hoverEnabled: true
        acceptedButtons: Qt.NoButton
        onPositionChanged: (mouse) => {
            const spots = root.checkWarningHotspots
            let hit = null
            for (let i = 0; i < spots.length; i++) {
                const s = spots[i]
                const dx = mouse.x - s.x, dy = mouse.y - s.y
                if (dx * dx + dy * dy <= s.r * s.r) { hit = s; break }
            }
            if (hit) {
                checkWarningTip.text = hit.text
                checkWarningTip.x = mouse.x + 8
                checkWarningTip.y = mouse.y + 8
                checkWarningTip.visible = true
            } else {
                checkWarningTip.visible = false
            }
        }
        onExited: checkWarningTip.visible = false
    }

    ToolTip {
        id: checkWarningTip
        visible: false
        timeout: -1
    }
}
